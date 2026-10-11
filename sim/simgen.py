"""物理一致性仿真器（正演）—— 先规定物理世界里发生了什么，再算 IMU 会读到什么。

为什么必须有它
--------------
真值目前只有用户肉眼估的 15 cm（本身 ±0.5 cm），拿它当基准分辨不出 3% 和 0.3%。
而"给加速度曲线加噪声再积分"是**逆运算造数据**：数据本身就来自算法要做的积分，
它验不出刚体几何（IMU 与测量头不重合）这类真正的建模错误 —— 那是当前最可疑的一项。

本仿真器方向相反：先定义测量头怎么走、手机怎么转，再正演出 IMU 读数。于是
  · 无噪声样例：任何正确的实现都必须精确复原真值 → 检验"实现正确性"；
  · 加噪声/零偏/比例因子/时差后跑大量随机种子 → 才谈得上"鲁棒性"。

坐标与符号（与 measure2.py / SlideMeasure.kt 完全一致，不能各说各话）
--------------------------------------------------------------------
  body  : 手机平放屏幕朝上时，屏幕右 +x，屏幕顶 +y，出屏幕朝天 +z
  world : 导航系，+z 朝天，滑动沿 +x
  加速度计输出比力 f_body = R^T·(a_I − g_vec)，g_vec = (0,0,−G)，故静止平放时 az = +G
  陀螺输出 ω_body（机体系角速度）
  q 表示 body→world：qrot(q, v_body) = v_world

刚体几何（本文件存在的核心理由）
--------------------------------
  IMU 在机体系原点，测量头在机体系固定偏移 r 处：
      p_H(t) = p_I(t) + R(t)·r        ⇒   p_I(t) = p_H(t) − R(t)·r
  IMU 的加速度不等于测量头的加速度；手机一转，即使测量头钉着不动，IMU 也会走一段圆弧。
"""
import math
import random

from pipeline import (qmul, qconj, qnorm, qrot, qrotT, q_from_euler,
                      vadd, vsub, vscale, vdot, vcross, vnorm, vunit)

G = 9.81


# ------------------------------------------------------------------ 配置

class SimCfg:
    """一台"仿真的手机"。数值默认取 K40 / lsm6dso 的实测量级（见 docs/calibration.md）。"""

    def __init__(self,
                 # ---- 时间 ----
                 hz=100.0,               # 名义采样率
                 dt_jitter_s=0.0,        # 采样间隔抖动 σ（真机实测 ~1 ms 量级）
                 gyr_offset_s=0.0,       # 陀螺相对加速度计的固定时差（正 = 陀螺晚）
                 drop_rate=0.0,          # 丢样概率
                 # ---- 加速度计 ----
                 acc_bias=(0.0, 0.0, 0.0),   # m/s²
                 acc_noise=0.010,            # m/s² 每样本白噪声 σ（K40 静止段实测 0.008~0.013）
                 acc_scale=(0.0, 0.0, 0.0),  # 比例因子误差（0.011 = +1.1%）
                 acc_misalign=0.0,           # 轴间非正交幅度（rad，~0.005 = 0.3°）
                 acc_bias_rw=0.0,            # m/s²/√s 零偏随机游走
                 acc_quant=0.0,              # 量化步长 m/s²（0.0024 ≈ ±8g/16bit）
                 acc_lp_hz=0.0,              # 抗混叠一阶低通截止（0 = 不滤波）
                 # ---- 陀螺 ----
                 gyr_bias=(0.0, 0.0, 0.0),   # rad/s
                 gyr_noise=0.002,            # rad/s 每样本 σ
                 gyr_scale=(0.0, 0.0, 0.0),
                 gyr_bias_rw=0.0,            # rad/s/√s（in-run instability）
                 # ---- 输出 ----
                 seed=0):
        self.__dict__.update(locals())
        del self.__dict__["self"]


def _misalign_matrix(cfg, rng):
    """I + 对角比例误差 + 对称/非对称非正交项。单位无量纲。"""
    sx, sy, sz = cfg.acc_scale
    m = rng
    e = cfg.acc_misalign
    return ((1 + sx, e * m(0), e * m(1)),
            (e * m(2), 1 + sy, e * m(3)),
            (e * m(4), e * m(5), 1 + sz))


def _apply_M(M, v):
    return (M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
            M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
            M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2])


# ------------------------------------------------------------------ 运动剖面

def smoothstep5(u):
    """0→1，一阶二阶导数在两端都为 0（所以起停都没有加速度冲击）。

    真机里手推的加速度不会在端点跳变，用这个比梯形速度剖面更贴近；
    副作用是峰值加速度更大（= 1.875·d/T² 的均值倍率），对噪声更友好但对零偏更苛刻，
    这本身就是"滑快一点"的物理理由。
    """
    u = min(1.0, max(0.0, u))
    return u * u * u * (10.0 - 15.0 * u + 6.0 * u * u)


def profile_s(t, t0, t1, d, pause=None):
    """测量头沿直线的标量位移 s(t)。pause=(pa,pb) 时中途停顿（两段各走一半）。"""
    if t <= t0:
        return 0.0
    if t >= t1:
        return d
    if pause is None:
        return d * smoothstep5((t - t0) / (t1 - t0))
    pa, pb = pause
    if t <= pa:
        return (d / 2) * smoothstep5((t - t0) / max(1e-9, pa - t0))
    if t >= pb:
        return (d / 2) + (d / 2) * smoothstep5((t - pb) / max(1e-9, t1 - pb))
    return d / 2


# ------------------------------------------------------------------ 姿态

def att_const(q0):
    return lambda t: q0


def att_slerp(q0, q1, t0, t1):
    """起止姿态不同（最考验姿态误差的一项）。"""
    def f(t):
        u = smoothstep5((t - t0) / max(1e-9, t1 - t0))
        return _slerp(q0, q1, u)
    return f


def _slerp(a, b, u):
    d = sum(a[i] * b[i] for i in range(4))
    if d < 0:
        b = tuple(-c for c in b)
        d = -d
    if d > 0.9995:
        return qnorm(tuple(a[i] + u * (b[i] - a[i]) for i in range(4)))
    th = math.acos(min(1.0, d))
    s = math.sin(th)
    wa = math.sin((1 - u) * th) / s
    wb = math.sin(u * th) / s
    return qnorm((wa * a[0] + wb * b[0], wa * a[1] + wb * b[1],
                  wa * a[2] + wb * b[2], wa * a[3] + wb * b[3]))


def att_wobble(q0, amp_deg=(1.0, 1.0, 1.0), hz=(1.1, 1.7, 0.7), phase=(0.0, 1.0, 2.0),
               drift_deg_per_s=(0.0, 0.0, 0.0), t0=0.0, env=None):
    """在基准姿态上叠手部抖动 + 缓慢零偏漂移。

    env(t)：幅度包络。真机里手抖只在**滑动期间**发生 —— 用户按下前后手机是放稳的
    （真机静止段 σ_a 只有 0.008~0.013 就是证据）。全程抖动是最坏情况，
    但拿它当唯一工况会把算法逼到"必须拒答"，反而看不出改进是否有效。
    所以两种都留：env=None 全程抖，env=滑动窗 只在滑动时抖。
    """
    a = tuple(math.radians(x) for x in amp_deg)
    w = tuple(2 * math.pi * x for x in hz)
    d = tuple(math.radians(x) for x in drift_deg_per_s)

    def f(t):
        tt = t - t0
        k = env(t) if env else 1.0
        r = k * (a[0] * math.sin(w[0] * tt + phase[0]) + d[0] * tt)
        p = k * (a[1] * math.sin(w[1] * tt + phase[1]) + d[1] * tt)
        y = k * (a[2] * math.sin(w[2] * tt + phase[2]) + d[2] * tt)
        return qnorm(qmul(q0, q_from_euler(r, p, y)))
    return f


def slide_envelope(tA, tB, edge=0.35):
    """滑动期间为 1，两侧 edge 秒内平滑降到 0。"""
    def f(t):
        if t <= tA - edge or t >= tB + edge:
            return 0.0
        if tA <= t <= tB:
            return 1.0
        if t < tA:
            return smoothstep5((t - (tA - edge)) / edge)
        return smoothstep5(((tB + edge) - t) / edge)
    return f


def compose_att(*fs):
    """把多个姿态项乘起来（右乘 = 在机体系里叠加）。"""
    def f(t):
        q = fs[0](t)
        for g in fs[1:]:
            q = qnorm(qmul(q, g(t)))
        return q
    return f


# ------------------------------------------------------------------ 正演

def imu_from_world(head_pos, att, r_body, t):
    """给定测量头位置函数、姿态函数、杠杆臂，返回该时刻 IMU 的比力与角速度。

    用中心差分求导数。步长取得比采样间隔小得多（1e-3 s vs 1e-2 s），
    截断误差 ~ h²·f''''/6 远小于我们要讨论的 cm 级误差；舍入误差 ~1e-16/h² 更小。
    """
    def p_I(s):
        return vsub(head_pos(s), qrot(att(s), r_body))

    h1, h2 = 1e-4, 2e-3
    # 角速度：ω_body = 2·vec(q⁻¹ ⊗ q̇)
    q = att(t)
    qp, qm = att(t + h1), att(t - h1)
    dq = tuple((qp[i] - qm[i]) / (2 * h1) for i in range(4))
    wq = qmul(qconj(q), dq)
    w_body = (2 * wq[1], 2 * wq[2], 2 * wq[3])

    # 加速度：p_I 的二阶中心差分
    p0, pp, pm = p_I(t), p_I(t + h2), p_I(t - h2)
    a_I = tuple((pp[i] - 2 * p0[i] + pm[i]) / (h2 * h2) for i in range(3))
    # 比力 f = a − g_vec = a + (0,0,G)
    f_nav = (a_I[0], a_I[1], a_I[2] + G)
    f_body = qrotT(q, f_nav)
    return f_body, w_body


def gen_timestamps(cfg, t_end, rng, offset_s=0.0):
    """名义均匀 + 抖动 + 丢样。返回 (ts_ns, t_s)。每条流独立生成 —— 真机就是两条独立流。"""
    dt0 = 1.0 / cfg.hz
    ts, t = [], 0.0
    while t <= t_end:
        if not (cfg.drop_rate and rng.random() < cfg.drop_rate):
            ts.append(int(round((t + offset_s) * 1e9)))
        t += dt0 + (rng.gauss(0, cfg.dt_jitter_s) if cfg.dt_jitter_s else 0.0)
    return ts, [x / 1e9 for x in ts]


def gen(cfg, head_pos, att, r_body, t_end, t0=0.0):
    """正演一次采集，返回与真机 CSV 同构的数据 + 真值。"""
    rng = random.Random(cfg.seed)
    Ma = _misalign_matrix(cfg, lambda k: rng.gauss(0, 1) * (1 if k < 6 else 1))
    Mg = ((1 + cfg.gyr_scale[0], 0, 0), (0, 1 + cfg.gyr_scale[1], 0),
          (0, 0, 1 + cfg.gyr_scale[2]))

    # 零偏随机游走：整段按 √t 漂移（in-run instability）
    def walk_per_sample(rate, dt):
        return rate * math.sqrt(dt) if rate else 0.0

    acc_ts, acc_t = gen_timestamps(cfg, t_end, rng, 0.0)
    gyr_ts, gyr_t = gen_timestamps(cfg, t_end, rng, cfg.gyr_offset_s)

    def sample(seq_t, M, bias, noise, rw, quant, lp_hz, is_acc):
        out = []
        ba = list(bias)
        prev = None
        dt_nom = 1.0 / cfg.hz
        for t in seq_t:
            if is_acc:
                f, w = imu_from_world(head_pos, att, r_body, t)
                true_v = f
            else:
                _, w = imu_from_world(head_pos, att, r_body, t)
                true_v = w
            ba = [ba[i] + rng.gauss(0, walk_per_sample(rw, dt_nom)) for i in range(3)]
            v = _apply_M(M, true_v)
            v = (v[0] + ba[0] + rng.gauss(0, noise),
                 v[1] + ba[1] + rng.gauss(0, noise),
                 v[2] + ba[2] + rng.gauss(0, noise))
            if lp_hz:  # 一阶低通（抗混叠 + 相位滞后，后者本身也是误差源）
                a = 1.0 / (1.0 + dt_nom * 2 * math.pi * lp_hz)
                v = tuple((a * v[i] + (1 - a) * (prev[i] if prev else v[i]))
                          for i in range(3))
            if quant:
                v = tuple(round(x / quant) * quant for x in v)
            prev = v
            out.append(v)
        return out

    acc = sample(acc_t, Ma, cfg.acc_bias, cfg.acc_noise, cfg.acc_bias_rw,
                 cfg.acc_quant, cfg.acc_lp_hz, True)
    gyr = sample(gyr_t, Mg, cfg.gyr_bias, cfg.gyr_noise, cfg.gyr_bias_rw,
                 0.0, 0.0, False)

    p0 = head_pos(t0)
    p1 = head_pos(t_end)
    return {
        "acc": acc, "gyr": gyr, "acc_ts": acc_ts, "gyr_ts": gyr_ts,
        "r_body": r_body,
        "headStart": p0, "headEnd": p1,
        "headNet": vnorm(vsub(p1, p0)),
        "headNetSigned": (p1[0] - p0[0]),
        "headPos": head_pos, "att": att,
        "tEnd": t_end, "t0": t0,
    }


def down_i_of(acc_ts, t_down):
    lo, hi = 0, len(acc_ts) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if acc_ts[mid] < t_down * 1e9:
            lo = mid + 1
        else:
            hi = mid
    return lo


# ------------------------------------------------------------------ 场景

def scenario(name, **kw):
    """返回一个 (cfg, head_pos, att, r_body, T, t0, truth) 的构造表。"""
    cfg = kw.pop("cfg", SimCfg())
    r_body = kw.pop("r_body", (0.0, 0.0, 0.0))
    T = kw.pop("T", 5.0)
    pre = kw.pop("pre", 1.5)      # 运动前静止
    post = kw.pop("post", 1.5)    # 运动后静止
    d = kw.pop("d", 0.15)
    t0 = 0.0
    tA, tB = pre, pre + kw.pop("move_s", 1.2)   # 运动起止
    t_end = tB + post

    q0 = qnorm(kw.pop("q0", (1.0, 0.0, 0.0, 0.0)))
    tilt = math.radians(kw.pop("tilt_deg", 0.0))
    if tilt:
        q0 = qnorm(qmul(q0, q_from_euler(tilt, 0.0, 0.0)))

    if name == "static":
        # 全程静止：手机不动。真值 0。专测"按着不动会不会造出距离"
        head = lambda t: (0.0, 0.0, 0.0)
        att = att_const(q0)
    elif name == "static_rotate":
        # 测量头钉死不动，手机绕它转（最狠的杠杆臂检验）：真值 0，IMU 走圆弧
        amp = kw.pop("rot_deg", 20.0)
        head = lambda t: (0.0, 0.0, 0.0)
        att = att_slerp(q0, qnorm(qmul(q0, q_from_euler(0.0, 0.0, math.radians(amp)))),
                        tA, tB)
    elif name == "jitter_only":
        # 只有手抖，测量头不移动。真值 0
        head = lambda t: (0.0, 0.0, 0.0)
        att = att_wobble(q0, kw.pop("amp_deg", (2, 2, 2)),
                         kw.pop("hz", (1.1, 1.7, 0.7)), t0=t0)
    elif name == "slide":
        # 标准滑动：测量头沿 +x 直线走 d 米，贴桌面（z 恒为 0）
        def head(t, _d=d, _a=tA, _b=tB):
            return (profile_s(t, _a, _b, _d, kw.get("pause")), 0.0, 0.0)
        att = att_const(q0)
    elif name == "slide_rotate":
        # 滑动 + 手机同时转动 —— 杠杆臂在真实场景里的样子
        def head(t, _d=d, _a=tA, _b=tB):
            return (profile_s(t, _a, _b, _d, kw.get("pause")), 0.0, 0.0)
        att = att_slerp(q0, qnorm(qmul(q0, q_from_euler(0.0, 0.0,
                        math.radians(kw.pop("rot_deg", 15.0))))), tA, tB)
    elif name == "slide_wobble":
        # 滑动 + 手抖（用户抱怨的"操作要求高"就是这一档）
        def head(t, _d=d, _a=tA, _b=tB):
            return (profile_s(t, _a, _b, _d, kw.get("pause")), 0.0, 0.0)
        env = slide_envelope(tA, tB) if kw.pop("only_during_slide", False) else None
        att = att_wobble(q0, kw.pop("amp_deg", (3, 3, 4)),
                         kw.pop("hz", (1.1, 1.7, 0.7)),
                         kw.pop("drift_deg_per_s", (0, 0, 2)), t0=t0, env=env)
    elif name == "slide_pivot":
        # 手机绕测量头倾斜（测量头贴着桌沿不动地转）之后再滑 —— 混合
        def head(t, _d=d, _a=tA, _b=tB):
            return (profile_s(t, _a, _b, _d, kw.get("pause")), 0.0, 0.0)
        att = compose_att(
            att_slerp(q0, qnorm(qmul(q0, q_from_euler(math.radians(
                kw.pop("roll_deg", 10.0)), 0.0, 0.0))), tA, tB),
            att_wobble((1.0, 0, 0, 0), kw.pop("amp_deg", (1, 1, 1)), t0=t0))
        att = lambda t: qnorm(qmul(att_slerp(q0, qnorm(qmul(q0, q_from_euler(
            math.radians(kw.get("roll_deg", 10.0)), 0.0, 0.0))), tA, tB)(t),
            att_wobble((1.0, 0, 0, 0), kw.get("amp_deg", (1, 1, 1)), t0=t0)(t)))
    else:
        raise ValueError(name)

    return dict(cfg=cfg, head_pos=head, att=att, r_body=r_body,
                t_end=t_end, t0=t0, tA=tA, tB=tB,
                truth=(0.0 if name in ("static", "static_rotate", "jitter_only") else d),
                name=name)


if __name__ == "__main__":
    # 冒烟：无噪声时，静止场景的 IMU 读数必须就是纯重力
    c = SimCfg(seed=1)
    s = scenario("static", cfg=c, r_body=(0.0, -0.09, 0.0), T=5.0)
    r = gen(s["cfg"], s["head_pos"], s["att"], s["r_body"], s["t_end"])
    a0 = r["acc"][10]
    print("acc[10] =", tuple(round(x, 6) for x in a0), " 期望 (0,0,9.81)")
    print("gyr[10] =", tuple(round(x, 6) for x in r["gyr"][10]))
    print("总样本 acc", len(r["acc"]), "gyr", len(r["gyr"]))

    s2 = scenario("static_rotate", cfg=c, r_body=(0.0, -0.09, 0.0), rot_deg=30)
    r2 = gen(s2["cfg"], s2["head_pos"], s2["att"], s2["r_body"], s2["t_end"])
    # 测量头没动，但 IMU 走了一段弧 —— 打印 IMU 的净位移量级
    print("测量头净位移（真值）: %.4f m" % r2["headNet"])
    f = [x for x in r2["acc"]]
    print("旋转中 acc 模长范围 %.3f ~ %.3f（静止时应恒为 9.81，转动时有离心/切向项）"
          % (min(vnorm(x) for x in f), max(vnorm(x) for x in f)))
