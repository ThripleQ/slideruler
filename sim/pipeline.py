"""slideruler 合成数据 + 算法管线仿真（纯标准库，无依赖）

目的：在真机之前，用已知真值的合成数据验证 docs/implementation-guide.md
阶段三~七的整条管线（姿态 → 去重力 → 一维投影 → ZUPT → 二次积分）
能不能把 0.8m 测进 5%~15%，以及这个结论对哪些参数敏感。

坐标系约定（沿用指南 §1.1）：
  body  : 手机平放时 屏幕右 = +x, 屏幕顶 = +y, 出屏幕朝天 = +z
  world : 滑动沿 +x, +z 指天
  加速度计读比力 f = a - g，g 为重力加速度矢量(0,0,-9.81) → 静止时 f_body ≈ (0,0,+9.81)
"""

import math
import random
import statistics

G = 9.81
DT = 0.01  # 100 Hz（第一版不含 jitter，见文末简化清单）

STATIC_PRE_S = 2.0
MOTION_S = 1.5
STATIC_POST_S = 1.2
D_TRUE = 0.8


# --------------------------------------------------------------------------
# 向量 / 四元数（q 表示 body -> world 的旋转）
# --------------------------------------------------------------------------
def vadd(a, b): return (a[0] + b[0], a[1] + b[1], a[2] + b[2])
def vsub(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def vscale(a, s): return (a[0] * s, a[1] * s, a[2] * s)
def vdot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def vcross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
def vnorm(a): return math.sqrt(vdot(a, a))
def vunit(a):
    n = vnorm(a)
    return vscale(a, 1.0 / n) if n > 1e-12 else (0.0, 0.0, 1.0)


def qmul(a, b):
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return (aw * bw - ax * bx - ay * by - az * bz,
            aw * bx + ax * bw + ay * bz - az * by,
            aw * by - ax * bz + ay * bw + az * bx,
            aw * bz + ax * by - ay * bx + az * bw)


def qnorm(q):
    n = math.sqrt(sum(c * c for c in q))
    return tuple(c / n for c in q) if n > 1e-12 else (1.0, 0.0, 0.0, 0.0)


def qconj(q): return (q[0], -q[1], -q[2], -q[3])


def qrot(q, v):
    """body -> world"""
    r = qmul(qmul(q, (0.0, v[0], v[1], v[2])), qconj(q))
    return (r[1], r[2], r[3])


def qrotT(q, v):
    """world -> body"""
    return qrot(qconj(q), v)


def q_from_euler(roll, pitch, yaw):
    cr, sr = math.cos(roll / 2), math.sin(roll / 2)
    cp, sp = math.cos(pitch / 2), math.sin(pitch / 2)
    cy, sy = math.cos(yaw / 2), math.sin(yaw / 2)
    return qnorm((cr * cp * cy + sr * sp * sy,
                  sr * cp * cy - cr * sp * sy,
                  cr * sp * cy + sr * cp * sy,
                  cr * cp * sy - sr * sp * cy))


def q_delta(w, dt):
    """角速度 w(rad/s) 持续 dt 的小旋转四元数"""
    return qnorm((1.0, w[0] * dt / 2.0, w[1] * dt / 2.0, w[2] * dt / 2.0))


def q_from_two_vectors(a_body, b_world):
    """求 q，使 R(q)·a_body == b_world（两向量均为单位向量）"""
    a, b = vunit(a_body), vunit(b_world)
    c = vcross(a, b)
    d = vdot(a, b)
    if d < -0.999999:  # 反向，绕任意垂直轴转 180°
        axis = vcross(a, (1.0, 0.0, 0.0))
        if vnorm(axis) < 1e-6:
            axis = vcross(a, (0.0, 1.0, 0.0))
        axis = vunit(axis)
        return (0.0, axis[0], axis[1], axis[2])
    return qnorm((1.0 + d, c[0], c[1], c[2]))


# --------------------------------------------------------------------------
# 配置
# --------------------------------------------------------------------------
class SensorCfg:
    """传感器与工况参数。默认值取自 docs/implementation-guide.md §1.3 与 §2。

    一台"设备" = 一组 SensorCfg。不同手机的差异就体现在这些字段上：
      acc_noise / acc_bias / acc_scale / gyro_* / dt
    """
    # 各字段的设备间典型差异（用于跨设备实验，见 run.py sec_device_matrix）
    def __init__(self,
                 acc_bias=0.05,        # m/s²，三轴独立常量零偏
                 acc_noise=0.02,       # m/s²，每样本白噪声 σ
                 gyro_bias=0.00873,    # rad/s（0.5 °/s），校准段能消掉的部分
                 gyro_noise=0.002,     # rad/s，每样本白噪声 σ
                 gyro_instab=3.5e-5,   # rad/s/√s（=0.002 °/s/√s），校准消不掉的慢漂
                 tilt_deg=2.0,         # 静态安装倾角
                 wobble_deg=1.0,       # 滑动中手部转动幅度
                 wobble_hz=1.2,
                 acc_scale=1.0,        # 加速度计标度因子（真机常在 0.98~1.02）
                 dt=DT,                # 真实采样间隔（秒）—— 真机实测 9.513 ms，别写死 10 ms
                 profile="sine"):      # "sine": 起始加速度最大(最友好)
                                       # "cycloid": 起始/终止加速度为 0(更像真实手推)
        self.acc_bias = acc_bias
        self.acc_noise = acc_noise
        self.gyro_bias = gyro_bias
        self.gyro_noise = gyro_noise
        self.gyro_instab = gyro_instab
        self.tilt_deg = tilt_deg
        self.wobble_deg = wobble_deg
        self.wobble_hz = wobble_hz
        self.acc_scale = acc_scale
        self.dt = dt
        self.profile = profile


class Params:
    """docs/implementation-guide.md §10 说的 PipelineParams：全部可调参数集中处。

    static_mode 是本次仿真新增的开关（指南原版只有 guide_magnitude 一种判据）：
      guide_magnitude : 指南 §7 判据一 |‖f‖-9.81| < 0.2  —— 本工况下失效，见 run.py
      linear_accel    : 改用去重力后的线性加速度模长（需要姿态，物理上正确）
      accel_variance  : 改用窗内加速度各轴标准差（不需要姿态和校准，最鲁棒）
      adaptive        : 阈值不写死，由本机静止段的实测噪声地板导出（自校准）
      oracle          : 直接用真值窗口（上界，用来隔离"检测器误差"与"算法误差"）

    assumed_dt 是"算法以为的采样间隔"。缺省 = 请求值 0.01 s；真机应按实测填入。
      位移 ∝ dt²，所以 dt 错 5% 距离就错 10% —— 这是最容易被忽略、代价最大的一项。

    use_scale_corr ⚠️ 有陷阱：它按 k̂ = ‖g_静态‖/9.81 反推标度。
      但单姿态的 ‖g_静态‖ 里混着"标度 k"和"零偏在重力方向的投影 b·ĝ"两项，无法分开
      （b_z 只要 0.05 m/s²，看起来就像 +0.5% 的标度）。因此**单姿态下不要开它**。
      只有当 g_mag 来自"多姿态标定"（≥4 个不共面朝向的球面拟合，见 fit_bias.py --selftest），
      或经独立手段已知 b≈0 时，才安全。
    """

    def __init__(self, alpha=0.05, static_acc_thr=0.2, static_gyr_thr=0.08,
                 static_win_s=0.3, use_attitude=True, use_zupt=True, use_pca=True,
                 reason_ok_rad=0.3, static_mode="guide_magnitude",
                 static_lin_thr=0.12, static_var_thr=0.05,
                 assumed_dt=DT, use_scale_corr=False,
                 adaptive_factor=1.6, adaptive_floor_lin=0.05, adaptive_floor_gyr=0.02,
                 cal_s=STATIC_PRE_S):
        self.alpha = alpha
        self.static_acc_thr = static_acc_thr
        self.static_gyr_thr = static_gyr_thr
        self.static_win_s = static_win_s
        self.use_attitude = use_attitude
        self.use_zupt = use_zupt
        self.use_pca = use_pca
        self.reason_ok_rad = reason_ok_rad
        self.static_mode = static_mode
        self.static_lin_thr = static_lin_thr
        self.static_var_thr = static_var_thr
        self.assumed_dt = assumed_dt
        self.use_scale_corr = use_scale_corr
        self.adaptive_factor = adaptive_factor
        self.adaptive_floor_lin = adaptive_floor_lin
        self.adaptive_floor_gyr = adaptive_floor_gyr
        self.cal_s = cal_s


# --------------------------------------------------------------------------
# 阶段〇：真值运动 + 传感器读数
# --------------------------------------------------------------------------
def gen_truth(cfg, seed=0, rng=None):
    rng = rng or random.Random(seed)
    dt = cfg.dt
    n_pre = int(round(STATIC_PRE_S / dt))
    n_mot = int(round(MOTION_S / dt))
    n_post = int(round(STATIC_POST_S / dt))
    N = n_pre + n_mot + n_post

    vpk = math.pi * D_TRUE / (2.0 * MOTION_S)
    w = math.pi / MOTION_S

    if cfg.profile == "cycloid":
        def _s(tau):
            return D_TRUE * (tau / MOTION_S - math.sin(2 * math.pi * tau / MOTION_S) / (2 * math.pi))
        def _a(tau):
            return (2.0 * D_TRUE / MOTION_S) * (math.pi / MOTION_S) * math.sin(2 * math.pi * tau / MOTION_S)
    else:  # "sine"：正弦速度剖面，起始加速度即峰值
        def _s(tau):
            return vpk / w * (1.0 - math.cos(w * tau))
        def _a(tau):
            return vpk * w * math.cos(w * tau)

    s = [0.0] * (N + 1)
    a_true = [0.0] * (N + 1)
    for i in range(N + 1):
        tau = (i - n_pre) * dt
        # 加速度要在滑动区间两端都取到真值采样，不能在边界处补 0
        if n_pre <= i <= n_pre + n_mot:
            a_true[i] = _a(tau)
        if i <= n_pre:
            s[i] = 0.0
        elif i < n_pre + n_mot:
            s[i] = _s(tau)
        else:
            s[i] = D_TRUE

    tilt = math.radians(cfg.tilt_deg)
    amp = math.radians(cfg.wobble_deg)
    ww = 2 * math.pi * cfg.wobble_hz

    q_true = []
    gyro_clean = []
    for i in range(N + 1):
        t = i * dt
        # 手部转动的包络必须连续（否则起点处姿态跳变，会伪造一个巨大的陀螺尖峰）
        u = (i - n_pre) * dt / MOTION_S
        env = math.sin(math.pi * u) ** 2 if 0.0 < u < 1.0 else 0.0
        k = amp * env
        roll = tilt + k * math.sin(ww * t)
        pitch = tilt * 0.6 + k * 0.7 * math.cos(ww * t * 0.8)
        q_true.append(q_from_euler(roll, pitch, 0.0))

    for i in range(N + 1):
        j = min(i + 1, N)
        dq = qmul(qconj(q_true[i]), q_true[j])
        gyro_clean.append((2.0 * dq[1] / dt, 2.0 * dq[2] / dt, 2.0 * dq[3] / dt))

    acc, gyr = [], []
    b = list(vunit(vscale(gyro_clean[0], 0.0)))  # 零偏随机游走的起点
    instability = [0.0, 0.0, 0.0]
    for i in range(N + 1):
        for a in range(3):
            instability[a] += cfg.gyro_instab * math.sqrt(dt) * rng.gauss(0, 1)
        f_world = (a_true[i], 0.0, G)  # 比力 = a - g, g=(0,0,-G)
        f_body = qrotT(q_true[i], f_world)
        # 传感器模型：acc_meas = k·f_body + b_acc + 白噪声（k 为标度因子）
        acc.append(tuple(cfg.acc_scale * f_body[a] + cfg.acc_bias + cfg.acc_noise * rng.gauss(0, 1)
                         for a in range(3)))
        gyr.append(tuple(gyro_clean[i][a] + cfg.gyro_bias + instability[a]
                         + cfg.gyro_noise * rng.gauss(0, 1) for a in range(3)))

    return {"N": N, "n_pre": n_pre, "n_mot": n_mot, "s": s, "a_true": a_true,
            "q_true": q_true, "acc": acc, "gyr": gyr, "dt": dt}


# --------------------------------------------------------------------------
# 阶段二：校准  —— 「算法校准」层
# --------------------------------------------------------------------------
def calibrate(truth, n_pre):
    """旧接口（debug/exp_* 脚本在用）：只返回陀螺零偏、重力参考单位向量、重力模长。"""
    acc, gyr = truth["acc"], truth["gyr"]
    win = range(0, n_pre)
    b_g = tuple(statistics.fmean(g[a] for g in [gyr[i] for i in win]) for a in range(3))
    a_mean = tuple(statistics.fmean(g[a] for g in [acc[i] for i in win]) for a in range(3))
    return b_g, vunit(a_mean), vnorm(a_mean)


class DeviceProfile:
    """从一次采集的「起始静止段」自动估计出的本机参数（免用户、每会话重算）。

    这就是「算法校准」的核心：不把任何阈值/零偏/重力常数写死，全部从本机当场数据导出。

      gyro_bias    ：静止段陀螺均值（校准段能消掉的那部分零偏）
      gravity_ref  ：静止段加速度均值方向（把 acc 常值零偏 + 安装倾角一并吸收）
      gravity_mag  ：静止段加速度均值模长 ≈ k·g_local（做去重力基准 & 标度修正）
      acc_sigma    ：静止段加速度噪声地板（用来导出静止检测阈值）
      gyro_sigma   ：静止段陀螺噪声地板

    为什么每会话重算、而不是出厂标定一次：陀螺零偏与噪声随温度漂移；安装倾角每次不同；
    而首尾静止段本来就强制存在 —— 校准数据是白送的，不花用户任何额外动作。
    """

    def __init__(self, dt, gyro_bias, gravity_ref, gravity_mag,
                 acc_sigma, gyro_sigma, n_static):
        self.dt = dt
        self.gyro_bias = gyro_bias
        self.gravity_ref = gravity_ref
        self.gravity_mag = gravity_mag
        self.acc_sigma = acc_sigma
        self.gyro_sigma = gyro_sigma
        self.n_static = n_static

    def scale_hat(self, g_nominal=G):
        """标度因子估计 k̂ = ‖g_meas‖ / g_nominal。

        ⚠️ 只用单姿态静止段时，这个 k̂ **不等于**加速度计标度 k：
           ‖g_meas‖ = ‖k·g_local·ĝ + b‖ ≈ k·g_local + b·ĝ，
           b 在重力方向的投影会直接伪装成标度（b=0.05 m/s² 就装成 +0.5%）。
        要真正分离 k 与 b，需要多姿态（≥4 个不共面朝向）球面拟合。
        另有不可消除的混淆：不知道当地 g 时，k̂ 里还混着当地纬度偏差（中国境内 ≤0.16%）。
        """
        return self.gravity_mag / g_nominal

    def summary(self):
        return (f"dt={self.dt*1000:.3f}ms   |g|={self.gravity_mag:.4f}   "
                f"k̂={self.scale_hat():.5f}   σ_a={self.acc_sigma:.4f}   "
                f"σ_g={self.gyro_sigma:.5f}   b_g={vnorm(self.gyro_bias):.5f} rad/s")


def build_profile(truth, p=None, g_nominal=G):
    """只用「起始静止段」估计本机参数。cal_s 秒，不依赖真值窗口。

    cal_s 默认 = 协议规定的静止时长（2.0 s）。真机上取保守一点更稳
    （用户可能没静满 2 s 就开推），但会略微增大零偏估计的随机误差。
    """
    dt = truth.get("dt", DT)
    cal_s = getattr(p, "cal_s", STATIC_PRE_S) if p is not None else STATIC_PRE_S
    n = max(4, int(round(cal_s / dt)))
    n = min(n, truth["N"])
    acc, gyr = truth["acc"], truth["gyr"]
    b_g = tuple(statistics.fmean(gyr[i][a] for i in range(n)) for a in range(3))
    a_mean = tuple(statistics.fmean(acc[i][a] for i in range(n)) for a in range(3))
    acc_sigma = max(statistics.pstdev([acc[i][a] for i in range(n)]) for a in range(3))
    gyro_sigma = max(statistics.pstdev([gyr[i][a] for i in range(n)]) for a in range(3))
    return DeviceProfile(dt, b_g, vunit(a_mean), vnorm(a_mean),
                         acc_sigma, gyro_sigma, n)


# --------------------------------------------------------------------------
# 阶段三：互补滤波姿态解算
# --------------------------------------------------------------------------
def solve_attitude(truth, b_g, g_ref_body, p, dt=None):
    acc, gyr, N = truth["acc"], truth["gyr"], truth["N"]
    dt = DT if dt is None else dt
    q0 = q_from_two_vectors(g_ref_body, (0.0, 0.0, 1.0))
    q = q0
    qs = [q]
    skipped = 0
    for i in range(N):
        w = vsub(gyr[i], b_g)
        q_pred = qnorm(qmul(q, q_delta(w, dt)))
        if not p.use_attitude:
            qs.append(q0)  # 对照：不做姿态跟踪，冻结在校准时刻
            continue
        g_pred_body = qrotT(q_pred, (0.0, 0.0, 1.0))
        e = vcross(g_pred_body, g_ref_body)
        if vnorm(e) > p.reason_ok_rad:      # 指南 §4 振动防护
            skipped += 1
            q = q_pred
        else:
            q = qnorm(qmul(q_pred, q_delta(vscale(e, -p.alpha), dt)))
        qs.append(q)
    return qs, skipped


def attitude_error_deg(qs, q_true):
    errs = []
    for q, qt in zip(qs, q_true):
        d = qnorm(qmul(qconj(qt), q))
        ang = 2.0 * math.acos(min(1.0, abs(d[0])))
        errs.append(math.degrees(ang))
    return errs


# --------------------------------------------------------------------------
# 阶段四：去重力
# --------------------------------------------------------------------------
def remove_gravity(truth, qs, g_mag=G):
    """去重力。g_mag 用 profile 实测的重力模长（而非写死 9.81）→ 顺带把标度误差摊平。"""
    a_lin = []
    for i in range(truth["N"] + 1):
        f_world = qrot(qs[i], truth["acc"][i])
        a_lin.append((f_world[0], f_world[1], f_world[2] - g_mag))
    return a_lin


# --------------------------------------------------------------------------
# 阶段五：PCA 一维投影
# --------------------------------------------------------------------------
def pca_axis(vectors):
    """幂迭代求协方差矩阵最大特征向量，返回 (u1, lambda1, lambda2)"""
    n = len(vectors)
    c = [[sum(v[a] * v[b] for v in vectors) / n for b in range(3)] for a in range(3)]
    u = (1.0, 0.0, 0.0)
    lam = 0.0
    for _ in range(200):
        nu = tuple(sum(c[a][b] * u[b] for b in range(3)) for a in range(3))
        lam = vnorm(nu)
        if lam < 1e-18:
            return (1.0, 0.0, 0.0), 0.0, 0.0
        u = vscale(nu, 1.0 / lam)
    # 收缩求第二特征值
    c2 = [[c[a][b] - lam * u[a] * u[b] for b in range(3)] for a in range(3)]
    u2 = (0.0, 1.0, 0.0)
    lam2 = 0.0
    for _ in range(200):
        nu = tuple(sum(c2[a][b] * u2[b] for b in range(3)) for a in range(3))
        lam2 = vnorm(nu)
        if lam2 < 1e-18:
            break
        u2 = vscale(nu, 1.0 / lam2)
    return u, lam, lam2


# --------------------------------------------------------------------------
# 阶段六：静止检测
# --------------------------------------------------------------------------
def detect_static(truth, b_g, p, a_lin=None, profile=None):
    acc, gyr, N = truth["acc"], truth["gyr"], truth["N"]
    dt = truth.get("dt", DT)
    wn = max(1, int(round(p.static_win_s / dt)))

    # adaptive 模式：阈值不写死，由本机「起始静止段」的实测噪声地板导出。
    thr_lin = p.static_lin_thr
    thr_gyr = p.static_gyr_thr
    if p.static_mode == "adaptive":
        if profile is None or a_lin is None:
            raise ValueError("adaptive 模式需要 profile 与 a_lin")
        nl = profile.n_static
        lin_floor = max(vnorm(a_lin[i]) for i in range(nl))
        gyr_floor = max(vnorm(vsub(gyr[i], b_g)) for i in range(nl))
        thr_lin = p.adaptive_factor * lin_floor + p.adaptive_floor_lin
        thr_gyr = p.adaptive_factor * gyr_floor + p.adaptive_floor_gyr

    flags = []
    for i in range(N + 1):
        lo = max(0, i - wn + 1)
        seg_g = [vnorm(vsub(gyr[k], b_g)) for k in range(lo, i + 1)]
        gyro_ok = max(seg_g) < thr_gyr

        if p.static_mode == "guide_magnitude":
            ok_a = max(abs(vnorm(acc[k]) - G) for k in range(lo, i + 1)) < p.static_acc_thr
        elif p.static_mode in ("linear_accel", "adaptive"):
            if a_lin is None:
                raise ValueError(f"{p.static_mode} 模式需要先算出去重力后的 a_lin")
            ok_a = max(vnorm(a_lin[k]) for k in range(lo, i + 1)) < thr_lin
        elif p.static_mode == "accel_variance":
            ok_a = max(statistics.pstdev([acc[k][a] for k in range(lo, i + 1)])
                       for a in range(3)) < p.static_var_thr
        elif p.static_mode == "oracle":
            ok_a = True
        else:
            raise ValueError(f"未知 static_mode: {p.static_mode}")

        flags.append(ok_a and gyro_ok)
    return flags


def motion_window(flags, truth=None, p=None):
    """第一个静止段结束 → 最后一个静止段开始"""
    n = len(flags)
    if p is not None and p.static_mode == "oracle":
        return truth["n_pre"], truth["n_pre"] + truth["n_mot"]
    start = next((i for i in range(n) if not flags[i]), None)
    end = next((i for i in range(n - 1, -1, -1) if not flags[i]), None)
    return start, end


# --------------------------------------------------------------------------
# 阶段七：ZUPT + 梯形二次积分
# --------------------------------------------------------------------------
def integrate_cond_zupt(a1, lo, hi, dt, bias_bound, sigma, k_sigma=3.0):
    """带「条件 ZUPT」的一维二次积分。

    原 ZUPT（`integrate_zupt`）无条件把窗口终点速度强制归零。这在手机"松手时还在动"
    时是救命的，但在手机已经停稳时是有害的：它把**残余零偏/低频漂移**攒出来的那点
    终点速度 v_end 当成真实速度扣掉，扣掉的位移 = |v_end|·T/2 —— 这个杠杆随窗口长度
    **线性**增长（而 v_end 本身又随 T 增长），所以长窗口下极易造成十几 % 的塌陷。

    实测（K40，20261011_090006，真值 15 cm，按住 4.40 s）：
        v_end = +0.0138 m/s（≈1.4 cm/s）
        无条件 ZUPT → 12.24 cm（−18.4%）；裸积分 → 15.27 cm（+1.8%）
    无条件 ZUPT 恰好扣掉 3.03 cm ≈ |v_end|·T/2 = 0.0138×4.40/2。

    判据：只有当 |v_end| 大到「用已知零偏量级解释不了」时才做去趋势。

        thr = bias_bound·T + k_sigma·σ_a·√(T·dt)

    第一项 = 已知残余零偏在 T 秒里最多能攒出的速度；第二项 = 白噪声在终点速度上的
    3σ 量级。|v_end| ≤ thr → 终点速度可用"停稳+零偏"解释 → 裸积分；
    否则 → 手机确实在动 → 去趋势。

    返回 (距离 m, v_end, 是否用了去趋势, thr)。
    """
    v_raw, p_raw = integrate(a1, lo, hi, dt=dt)
    T = (hi - lo) * dt
    v_end = v_raw[-1]
    thr = bias_bound * T + k_sigma * sigma * math.sqrt(max(T * dt, 0.0))
    used = abs(v_end) > thr
    if not used:
        return abs(p_raw[-1]), v_end, False, thr
    da = v_end / T if T > 1e-9 else 0.0
    a1c = [a1[i] - da if lo <= i <= hi else a1[i] for i in range(len(a1))]
    _, p_z = integrate(a1c, lo, hi, dt=dt)
    return abs(p_z[-1]), v_end, True, thr


def integrate(acc1, lo, hi, dt=DT):
    v = [0.0] * (hi - lo + 1)
    p = [0.0] * (hi - lo + 1)
    for k in range(1, hi - lo + 1):
        i = lo + k
        v[k] = v[k - 1] + (acc1[i - 1] + acc1[i]) / 2.0 * dt
        p[k] = p[k - 1] + (v[k - 1] + v[k]) / 2.0 * dt
    return v, p


def run_pipeline(truth, cfg, p, seed=0):
    # ---- 阶段二：自校准（只用起始静止段，免用户）----
    prof = build_profile(truth, p)
    dt_used = p.assumed_dt            # 算法"以为"的采样间隔（真机应填实测值）
    b_g, g_ref_body = prof.gyro_bias, prof.gravity_ref

    qs, skipped = solve_attitude(truth, b_g, g_ref_body, p, dt=dt_used)

    g_used = prof.gravity_mag if p.use_scale_corr else G
    a_lin = remove_gravity(truth, qs, g_mag=g_used)
    if p.use_scale_corr:
        k = prof.scale_hat()
        a_lin = [vscale(a_lin[i], 1.0 / k) for i in range(len(a_lin))]

    flags = detect_static(truth, b_g, p, a_lin=a_lin, profile=prof)
    lo, hi = motion_window(flags, truth, p)
    if lo is None or hi is None or hi - lo < 10:
        return {"ok": False, "reason": "静止检测失败，未找到滑动段"}

    vecs = a_lin[lo:hi + 1]
    u, lam1, lam2 = pca_axis(vecs)
    axis = u if p.use_pca else (1.0, 0.0, 0.0)
    # 主轴取正向（朝滑动方向）
    if vdot(u, (1.0, 0.0, 0.0)) < 0:
        axis = vscale(axis, -1.0) if p.use_pca else axis
    a1 = [vdot(a_lin[i], axis) for i in range(len(a_lin))]

    v_raw, p_raw = integrate(a1, lo, hi, dt=dt_used)
    dist_raw = abs(p_raw[-1])

    if not p.use_zupt:
        return {"ok": True, "dist": dist_raw, "dist_raw": dist_raw, "v_err": v_raw[-1],
                "skipped": skipped, "lam_ratio": (lam1 / lam2) if lam2 > 0 else float("inf"),
                "axis": axis, "lo": lo, "hi": hi, "dist_zupt": None,
                "dt": dt_used, "profile": prof, "g_mag": prof.gravity_mag}

    # 线性漂移扣除：把终点残余速度反推成恒定假加速度，全程扣掉重积一次
    T_mot = (hi - lo) * dt_used
    da = v_raw[-1] / T_mot
    a1_corr = [a1[i] - da if lo <= i <= hi else a1[i] for i in range(len(a1))]
    v_z, p_z = integrate(a1_corr, lo, hi, dt=dt_used)
    dist_zupt = abs(p_z[-1])

    att_err = attitude_error_deg(qs, truth["q_true"])
    return {"ok": True, "dist": dist_zupt, "dist_zupt": dist_zupt, "dist_raw": dist_raw,
            "v_err": v_raw[-1], "skipped": skipped,
            "lam_ratio": (lam1 / lam2) if lam2 > 0 else float("inf"),
            "axis": axis, "lo": lo, "hi": hi,
            "att_err_static": statistics.fmean(att_err[:truth["n_pre"]]),
            "att_err_motion": statistics.fmean(att_err[lo:hi + 1]),
            "g_mag": prof.gravity_mag, "q_est": qs, "a_lin": a_lin,
            "v_raw": v_raw, "p_raw": p_raw, "p_zupt": p_z, "a1": a1,
            "dt": dt_used, "profile": prof}
