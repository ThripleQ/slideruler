"""测距链路 v2 —— 先把「静止」和「运动」分开，再在两个静止锚点之间积分。

为什么推倒重来
--------------
v1 用的是「条件 ZUPT」：|v_end| ≤ 阈值时走裸积分、否则去趋势。这个判据**自相矛盾** ——
它拿零偏估计 b 去判断"终点速度是不是零偏造成的"，而零偏正是误差的来源。于是任何一个
静止窗口都必然通过检验（|v_end| = b·T 恒等于阈值里的第一项），必然输出 ½·b·T² 的假位移。

实测：`20261011_100631`，手机 3.64 s 全程静止（窗口内静止样本 100%），
|b|=0.035 m/s² → 裸积分 **18.20 cm**。这正是用户报的"按着不动也测出很长一段距离"。

v2 的三条原则
-------------
1. **静止检测只用原始数据**（加速度局部方差 + 陀螺模长），不依赖姿态、不依赖零偏。
   这样它不会和它要保护的东西互相污染。
2. **窗口 = 真实运动段**，由前后两段验证过的静止夹出来 —— 不是用户按住的时长。
   静止段里的真实位移恒为 0，却贡献完整的零偏误差（∝T²），所以多算一秒纯亏。
   用户"按住→停顿→推→停稳→松手"里的停顿和停稳，本该被切掉。
3. **两端都验证过静止 ⇒ 无条件去趋势**，而且把两侧静止段里的**所有**样本一起当锚点
   （双锚点 ZUPT），而不是只钉两个端点。噪声按 √N 平均掉。

坐标与量纲：acc 单位 m/s²，gyr 单位 rad/s，dt 单位 s，输出米。
"""
import math
import os
import statistics

from pipeline import (G, q_delta, q_from_two_vectors, qnorm, qrot, qrotT, qmul,
                       vadd, vcross, vnorm, vscale, vsub, vunit)

# ------------------------------------------------------------------ 参数

STATIC_WIN_S = 0.20      # 静止检测的滑动窗长（边界要锐）
SLOW_WIN_S = 0.80        # 慢动作检测的滑动窗长（白噪声的 σ 不随窗长变，慢信号却会涨）
STATIC_FACTOR = 4.0      # 静止阈值 = 本机噪声地板 × 该系数
STATIC_FLOOR = 0.010     # 噪声地板的硬下限（m/s²），防止极端记录把阈值压到噪声里
STATIC_MIN_S = 0.20      # 一段静止至少这么久才敢当锚点（约束的是整段静止，不是裁进窗口的余量）
GYRO_FLOOR = 0.010       # 陀螺静止阈值的硬下限（rad/s）
CAL_MIN_S = 0.25         # 校准段最短时长（0.25 s 足够：陀螺零偏估计误差 σ/√N ≈ 1.4e-4 rad/s）
CAL_MAX_S = 2.0          # 校准段最长（只用运动前这一段的末尾，别拿太旧的数据）
CAL_STATIC_FRAC = 0.95   # 校准段里静止样本的最低占比（混着运动样本的均值既不是
#                          零偏也不是运动加速度，拿它当重力方向会漏一大块。实测 100339：
#                          |g|=9.43 应为 9.91、σ=2.01，v1 据此报了 168.62 cm）
MARGIN_S = 0.15          # 运动段两侧各留的静止余量（落在静止段内，所以无害）
TAIL_S = 2.5             # UP 之后最多再看这么久，用来找运动结束后的静止段
PRE_S = 2.0              # DOWN 之前最多往回看这么久，用来找校准段
GAP_S = 0.15             # 两段运动之间隔这么久以内的静止，算同一段运动
# 运动区间的选法。见 measure(merge_mode=...) 的注释。
# "press" = 按钮按住期间发生过的运动取并集（v3 默认，不再靠静止检测切段）
# "overlap" = 取与按钮窗重叠最多的那一段（v2 行为，会切碎滑动、且惩罚滑过松手时刻的段）
MERGE_MODE = "press"
ALPHA = 0.05             # 互补滤波系数
GATE_RAD = 0.30          # 加速度修正的门限（预测/实测重力方向夹角）—— 仅 g_mag 未知时用
GATE_MAG_FRAC = 0.25     # 加速度修正的门限：|acc| 偏离 |g| 超过这个比例就不修正。
#   为什么换成模长检验：方向门限会**自锁**。20261011_143734 实测，姿态误差从
#   一开始就有 19.5°（|e|=0.335 > 0.30），于是全程"跳过修正"，误差再也拉不回来。
#   门限本意是"别拿运动加速度当重力"，但静止与否已由 flags 判定；
#   模长检验更能表达同一件事：静止时 |acc| 就该等于 |g|，与方向无关。
PLATEAU_GRID = 5         # 平台宽度扫描的点数
QUAD = False             # 残差加速度是否允许线性漂移（见 zupt_two_anchor 的说明）
# 「修正量相对答案」的上限：|裸积分 − 去趋势| 超过 max(答案, 该下限) 就报"模型敏感"。
# 0.02 m 这个下限是为了避免小距离时被绝对噪声误报。
MODEL_GAP_FLOOR_M = 0.02


# ------------------------------------------------------------------ 小工具

def rolling_max_std(rows, dt, win_s=STATIC_WIN_S):
    """逐样本局部标准差（居中窗，逐轴取最大）。居中 → 运动边界不滞后。"""
    n = len(rows)
    h = max(1, int(win_s / dt) // 2)
    out = [0.0] * n
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        c = b - a
        best = 0.0
        for x in range(3):
            m = 0.0
            for k in range(a, b):
                m += rows[k][x]
            m /= c
            s = 0.0
            for k in range(a, b):
                d = rows[k][x] - m
                s += d * d
            best = max(best, math.sqrt(s / c))
        out[i] = best
    return out


def rolling_max_norm(vecs, dt, win_s=STATIC_WIN_S):
    n = len(vecs)
    h = max(1, int(win_s / dt) // 2)
    out = [0.0] * n
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        out[i] = max(vnorm(vecs[k]) for k in range(a, b))
    return out


def percentile(xs, q):
    s = sorted(xs)
    if not s:
        return 0.0
    return s[min(len(s) - 1, max(0, int(q * len(s))))]


def static_flags(acc, gyr, gb, dt, floor_a=None, floor_g=None):
    """静止标记 + 用到的阈值（全部来源可追溯）。

    floor_a / floor_g 给了就不用 5% 分位数自举。自举有个致命前提：
    数据里得真有 ≥5% 的样本是静止的。手机全程手抖时这个前提不成立 ——
    那 5% 是抖得最轻的时刻，阈值被抬到抖动幅度之上 4 倍，
    连真正的滑动都判成静止（仿真 slide_wobble：15 cm 被报成 0.20 cm）。
    器件噪声地板是可以离线标定一次的常数，比每次自举可靠。
    """
    sd = rolling_max_std(acc, dt)
    floor = percentile(sd, 0.05) if floor_a is None else floor_a
    thr_a = max(STATIC_FACTOR * floor, STATIC_FLOOR)
    om = [vsub(gyr[i], gb) for i in range(len(gyr))]
    on = rolling_max_norm(om, dt)
    ofloor = percentile(on, 0.05) if floor_g is None else floor_g
    thr_g = max(STATIC_FACTOR * ofloor, GYRO_FLOOR)
    flags = [sd[i] <= thr_a and on[i] <= thr_g for i in range(len(acc))]
    return flags, floor, thr_a, ofloor, thr_g


def runs_of(flags, want, min_len):
    """返回 [start, end) 形式的连续段。"""
    out = []
    i = 0
    n = len(flags)
    while i < n:
        if flags[i] == want:
            j = i
            while j + 1 < n and flags[j + 1] == want:
                j += 1
            if j - i + 1 >= min_len:
                out.append((i, j + 1))
            i = j + 1
        else:
            i += 1
    return out


# ------------------------------------------------------------------ 姿态

def solve_attitude_v2(acc, gyr, gb, g_ref_body, flags, dt, alpha=ALPHA, gate=GATE_RAD,
                      use_flags=None, i0=0, g_mag=None):
    """互补滤波。**加速度修正在静止段才开**，运动中纯陀螺外推。

    i0    从第 i0 个样本才开始积分（之前的填初始姿态）。
          必须传校准段起点：初始姿态取自校准段的重力方向，若从样本 0 积起，
          而手机在校准段之前转过（实测 143734 转了 21°），这个初始姿态对
          样本 0 就是错的，而且会一路带进窗口 —— 静止段残余达 3.3 m/s²。
    g_mag 给了就用模长检验代替方向门限（见 GATE_MAG_FRAC 的说明）。

    use_flags：姿态修正的门控（默认 = flags）。两者可以不同 —— 这正是 v4 的关键改动：
      flags     问的是「这段能不能当 ZUPT 锚点」（要保守，宁缺毋滥）
      use_flags 问的是「此刻的 acc 方向能不能当重力方向用」（要严格：只要手机在动，
                 acc = 重力 + 运动加速度 + 杠杆臂的切向/向心项，它的方向**不是**重力方向）
    手抖 3° 经过 9 cm 杠杆臂就有 0.22 m/s² 的额外加速度，相当于把重力方向拧了 1.3°，
    再拿它去"修正"姿态，等于把误差写进姿态里 —— 重力泄漏 0.22 m/s²，
    在 2 s 窗口上积分两步就是 10 cm 量级。抖动期间**必须**纯陀螺外推。
    运动段由按钮窗界定（外部事实，不受数据自举污染）。

    v1 的参考是 `g_ref_body`（校准时刻的体重力方向，体坐标里固定不动）——
    等于让滤波器主动抵抗手机的真实旋转：持续角速率 ω 下稳态姿态误差 ≈ ω·dt/α
    （α=0.05、dt=9.5 ms 时 ≈ 0.19·ω），ω=5°/s 就有 1° 误差、重力泄漏 0.17 m/s²。
    v2 换成当前**实测**加速度方向（标准 Mahony），手机转它就跟着转，没有这个稳态误差。
    """
    ok = flags if use_flags is None else use_flags
    q = q_from_two_vectors(g_ref_body, (0.0, 0.0, 1.0))
    # qs[i] 必须对应 acc[i]，所以前 i0+1 个都是初始姿态（acc[i0] 时刻的姿态就是它）
    i0 = max(0, min(i0, max(0, len(acc) - 1)))
    qs = [q] * (i0 + 1)
    for i in range(i0, len(acc) - 1):
        q_pred = qnorm(qmul(q, q_delta(vsub(gyr[i], gb), dt)))
        if ok[i]:
            e = vcross(qrotT(q_pred, (0.0, 0.0, 1.0)), vunit(acc[i]))
            good = (abs(vnorm(acc[i]) - g_mag) <= GATE_MAG_FRAC * g_mag
                    if g_mag else vnorm(e) <= gate)
            q = (q_pred if not good
                 else qnorm(qmul(q_pred, q_delta(vscale(e, -alpha), dt))))
        else:
            q = q_pred
        qs.append(q)
    while len(qs) < len(acc):
        qs.append(qs[-1] if qs else q)
    return qs


def to_linear(acc, qs, g_mag):
    out = []
    for i in range(len(acc)):
        f = qrot(qs[i], acc[i])
        out.append((f[0], f[1], f[2] - g_mag))
    return out


# ------------------------------------------------------------------ 积分

def integrate_from(a_lin, i0, i1, dt):
    """从 i0 积到 i1（含），返回 (v, p)，下标 0..i1-i0。"""
    m = i1 - i0
    v = [[0.0, 0.0, 0.0] for _ in range(m + 1)]
    p = [[0.0, 0.0, 0.0] for _ in range(m + 1)]
    for k in range(1, m + 1):
        i = i0 + k
        for x in range(3):
            v[k][x] = v[k - 1][x] + 0.5 * (a_lin[i - 1][x] + a_lin[i][x]) * dt
            p[k][x] = p[k - 1][x] + 0.5 * (v[k - 1][x] + v[k][x]) * dt
    return v, p


def zupt_two_anchor(a_lin, w0, w1, run0, run1, dt, quad=False):
    """双锚点 ZUPT：把两侧静止段内**所有**样本的速度均值拉到 0。

    设 v(t) = v_true(t) + e(t)，e 是残差加速度的积分。A、B 由两个静止段的
    速度均值解出（取段内平均时刻/平均 t²、平均速度）。比只钉两个端点稳 √N。

    quad=False：e(t) = A·t + B            —— 2 个参数/轴，假设残差加速度是常量
    quad=True ：e(t) = A·t + B·t²/2       —— 2 个参数/轴，允许残差加速度线性漂移
      两段静止各给 3 个约束（速度矢量=0），正好够定 6 个未知数。
      代价：窗口短时矩阵接近奇异，噪声会被放大 —— 只在长窗口才划得来。
    """
    v, _ = integrate_from(a_lin, w0, w1, dt)

    # 锚点段裁到窗口内（平台扫描会传入比窗口更宽的段）
    run0 = (max(run0[0], w0), min(run0[1], w1 + 1))
    run1 = (max(run1[0], w0), min(run1[1], w1 + 1))
    if run0[1] <= run0[0]:
        run0 = (w0, w0 + 1)
    if run1[1] <= run1[0]:
        run1 = (w1, w1 + 1)

    def run_moments(v0, v1):
        idx = [i - w0 for i in range(v0, v1)]
        c = len(idx)
        t = sum(i * dt for i in idx) / c
        t2 = sum((i * dt) ** 2 for i in idx) / c
        vv = [sum(v[i][x] for i in idx) / c for x in range(3)]
        return t, t2, vv

    t0, t02, V0 = run_moments(*run0)
    t1, t12, V1 = run_moments(*run1)
    AB = []
    for x in range(3):
        a11, a12, a21, a22 = (t0, t02 / 2.0, t1, t12 / 2.0) if quad else \
                             (t0, 1.0, t1, 1.0)
        det = a11 * a22 - a12 * a21
        if abs(det) < 1e-12:
            AB.append((0.0, 0.0))
        else:
            AB.append(((V0[x] * a22 - a12 * V1[x]) / det,
                       (a11 * V1[x] - V0[x] * a21) / det))

    def corr(k):
        t = k * dt
        if quad:
            return [AB[x][0] * t + AB[x][1] * t * t / 2.0 for x in range(3)]
        return [AB[x][0] * t + AB[x][1] for x in range(3)]

    vc = [[0.0, 0.0, 0.0] for _ in v]
    p = [[0.0, 0.0, 0.0] for _ in v]
    for k in range(len(v)):
        c = corr(k)
        for x in range(3):
            vc[k][x] = v[k][x] - c[x]
            if k:
                p[k][x] = p[k - 1][x] + 0.5 * (vc[k - 1][x] + vc[k][x]) * dt
    dist = vnorm(p[-1])
    # 锚点漂移：**修正前**速度在两侧静止段里的均值模长 —— 也就是"要去掉多少速度"。
    # （早先这里算的是修正**后**的均值，那是构造上恒等于 0 的死字段。）
    # 它和 |b|·T 同量级时才说明残差是常量；远大于 |b|·T 说明静止段里其实还在动。
    anchor_drift = max(vnorm(V0), vnorm(V1))
    # 段内速度的散布（理想：静止段里速度应当平，散的是噪声）
    spread0 = statistics.pstdev([vnorm(vc[i - w0]) for i in range(*run0)])
    return dist, vc, p, anchor_drift, max(spread0, statistics.pstdev(
        [vnorm(vc[i - w0]) for i in range(*run1)]))


def anchor_gyro_tremor(gyr, run, dt):
    """静止段里的**角抖动**（rad/s）—— 它不是陀螺零偏不稳定性，别再那么叫。

    做法：把静止段切成 4 段，相邻段均值之差的模长 /√2，就是该时长尺度上
    "角速度在这段里变了多少"。切 4 段而不是更多，是为了让 τ 尽量长：
    白噪声在段均值里按 1/√h 衰减，段越短估计越大。

    ## 它到底量的是什么（2026-10-16，用库里那段 76 s 静止录制查清）

    `20261011_085650` 是一段 76 s 的纯静止录制（meta 里按住时长 0 ms，没做滑动）。
    在那一段上：
      · 5 s 窗口的零偏均值全部落在 ±3e-5 rad/s 内，随机散布、无单调漂移；
      · Allan deviation 从 10 ms 到 9.7 s **一路单调下降**（0.00052 → 7.8e-6 rad/s），
        全程没有平台（flicker 的斜率 0）也没有回升（速率随机游走的 +0.5）。
    ⟹ **这台 K40 的零偏不游走。** 该静止段上本函数只给 1.4e-5 rad/s。

    可它在别的采集上给的是 0.0001~0.17 —— 跨 3 个数量级，而且**按场次干净分层**：
      桌面受控场次 100444/100459/100529 → 0.00012~0.00041（≈ 器件噪声底）
      "姿态随意"手持场次 143719/143725/143734 → 0.0123~0.0200（差 30~170 倍）
    ⟹ 它量的是**起手前那一下手抖不抖**，即端点参考稳不稳，不是器件漂移。

    `attitude_drift_budget`（(1/6)·g·σ·√(T/τ)·T³）已删，理由两条：
      ① 物理前提被上面那段数据证伪（零偏既不游走、σ 也不随 τ 增长）；
      ② 它算出的 T³ 预算与实际误差**反相关** —— 预算 5.4 cm 的 143734 错得
         最狠（−90%），预算 35.7 cm 的 143725 反而 +5%。一个被自己的数据
         违反 2.5 倍的"上界"不是上界。

    返回 (σ rad/s, τ s)。τ 必须一起报 —— 同一个角抖动在 τ=76 ms 与 τ=304 ms
    上不是同一个数，脱离 τ 谈它是没意义的。
    """
    a, b = run
    n = b - a
    if n < 16:
        return None, None
    h = max(4, n // 4)
    k = n // h
    if k < 2:
        return None, None
    segs = [tuple(statistics.fmean(gyr[i][x] for i in range(a + j * h, a + (j + 1) * h))
                  for x in range(3)) for j in range(k)]
    d = [vnorm(vsub(segs[j], segs[j + 1])) for j in range(k - 1)]
    return statistics.fmean(d) / math.sqrt(2.0), h * dt


def principal_dir(pts, iters=60):
    """点云的第一主成分方向（幂法）。3×3 的协方差矩阵不需要通用特征分解。

    用它而不是"端点连线"，是因为端点只有两个样本，噪声直接进方向；
    主方向用的是整条轨迹，端点噪声被 √N 摊薄。
    """
    c = len(pts)
    if c < 2:
        return (1.0, 0.0, 0.0)
    mu = tuple(sum(p[i] for p in pts) / c for i in range(3))
    C = [[0.0] * 3 for _ in range(3)]
    for p in pts:
        d = (p[0] - mu[0], p[1] - mu[1], p[2] - mu[2])
        for i in range(3):
            for j in range(3):
                C[i][j] += d[i] * d[j]
    v = (1.0, 1.0, 1.0)
    for _ in range(iters):
        w = tuple(sum(C[i][j] * v[j] for j in range(3)) for i in range(3))
        nv = vnorm(w)
        if nv < 1e-18:
            return vunit(v)
        v = tuple(x / nv for x in w)
    return vunit(v)


# ------------------------------------------------------------------ 主流程

def measure(acc, gyr, down_i, up_i, dt=None, tail_s=TAIL_S, verbose=False, quad=QUAD,
            merge_mode=MERGE_MODE, head_r=None,
            floor_a=None, floor_g=None, g_mag_fix=None, gate_motion=True):
    """merge_mode 决定「哪一段运动算这次测量」：
      "overlap" 取与按钮窗重叠最多的那一组（v2 原行为，见下面注释里的缺陷）
      "press"   取按钮按住期间发生过的所有运动的**并集**（v3，默认）
                —— 不再靠静止检测去切运动段，所以匀速段不会被当成"停了"把滑动切碎。
      "button"  **直接以按钮窗口为运动区间**（v4）。press 仍要靠静止检测认出"哪里在动"，
                而静止检测的噪声地板是从本次数据自举出来的（5% 分位数）——
                手机全程手抖时，那 5% 就是抖得最轻的时刻，阈值被抬到抖动幅度之上，
                连真正的滑动都判成"静止"。仿真 slide_wobble 里运动段只剩 0.07 s、
                15 cm 报成 0.20 cm，就是这么来的。按钮时刻是用户给的外部事实，
                不受数据自举污染，所以边界应该由它定，静止检测退回只管找锚点。

    floor_a / floor_g：器件噪声地板（m/s²、rad/s）。给了就用它、不再自举。
    """
    n = min(len(acc), len(gyr))
    if dt is None:
        dt = 0.01
    res = {"ok": False, "reason": None, "n": n}

    # ---- 1. 静止检测（只用原始数据）----
    #       陀螺零偏先用「最安静的 0.5 s」粗估 —— 只为了把 ω 的直流分量去掉，精度要求很低。
    pre = max(20, min(n // 4, int(0.5 / dt)))
    sd0 = rolling_max_std(acc, dt)
    i_quiet = min(range(n), key=lambda i: sd0[i])
    qa = max(0, i_quiet - pre // 2)
    qb = min(n, i_quiet + pre // 2 + 1)
    gb = tuple(statistics.fmean(gyr[i][a] for i in range(qa, qb)) for a in range(3))
    flags, floor, thr_a, ofloor, thr_g = static_flags(acc, gyr, gb, dt, floor_a, floor_g)
    res.update(floor=floor, thr_a=thr_a, ofloor=ofloor, thr_g=thr_g, quiet=(qa, qb))

    min_run = max(3, int(STATIC_MIN_S / dt))
    static_runs = runs_of(flags, True, min_run)

    # ---- 2. 搜索区：DOWN 前 2 s → UP 后 tail_s ----
    # 按钮时段必须有正长度。DOWN==UP 时所有运动组与按钮窗的重叠都是 0，
    # 下面的"取重叠最大的一组"就退化成"取时间最早的一组" —— 会去测**按住之前**
    # 那一小段（用户挪手机）而不是这次滑动。实测 100444 上它会报 0.06 cm 而不是 15.4 cm：
    # 数字看着无害，但它答的是另一个问题。宁可拒答。
    if up_i <= down_i:
        return {**res, "reason": "按下与松手时刻异常（松手不晚于按下）"}
    reg_lo = max(0, down_i - int(PRE_S / dt))
    reg_hi = min(n - 1, up_i + int(tail_s / dt))
    gap = max(2, int(GAP_S / dt))

    mot = [i for i in range(reg_lo, reg_hi + 1) if not flags[i]]
    if not mot:
        # 0.2 s 的局部标准差看不出一场很慢的滑动（15 cm 滑 4 s，峰值加速度只有
        # 0.059 m/s²，0.2 s 窗内的变化量 ≈0.006 —— 埋在噪声里）。但它在 0.8 s 的
        # 尺度上看得见。这里不猜距离，只把"其实移动了"这件事说出来，让用户滑快一点。
        sd_slow = rolling_max_std(acc, dt, SLOW_WIN_S)
        slow_peak = max(sd_slow[i] for i in range(reg_lo, reg_hi + 1))
        slow_thr = max(STATIC_FACTOR * percentile(sd_slow, 0.05), STATIC_FLOOR)
        if slow_peak > slow_thr:
            return {**res, "slowMotion": True, "slowPeak": slow_peak, "slowThr": slow_thr,
                    "reason": f"移动太慢，检测器分辨不出（{SLOW_WIN_S}s 尺度上只有 "
                              f"{slow_peak:.4f} / 阈值 {slow_thr:.4f} m/s²）—— 请快速滑完（1~2 秒）"}
        return {**res, "ok": True, "distanceM": 0.0, "noMotion": True,
                "imuDistM": 0.0, "axisDistM": 0.0, "leverM": 0.0,
                "straightness": 0.0, "straightRms": 0.0,
                "reason": "窗口内没有任何运动", "plateau": (0.0, 0.0),
                "slowPeak": slow_peak, "slowThr": slow_thr,
                "staticRuns": static_runs, "endHitBound": False}

    # 按间隙分组
    groups = [[mot[0], mot[0]]]
    for i in mot[1:]:
        if i - groups[-1][1] <= gap:
            groups[-1][1] = i
        else:
            groups.append([i, i])

    def overlap(g):
        return max(0, min(g[1], up_i) - max(g[0], down_i))

    if merge_mode == "press":
        # ---- v3：按钮按住期间发生过的运动，全部算同一次测量，取并集 ----
        # 不能靠静止检测去切运动段。**匀速运动在加速度计上和静止无法区分**
        # （伽利略相对性原理：加速度计测的是比力，匀速时比力就是重力），
        # 所以一次平稳滑动的中间段必然被判成"静止"，一次完整滑动被切成好几段。
        # 真机证据 122416：一次按住里被切成 5 段（0.20/0.29/0.14/0.26/1.19 s），
        # 只取其中一段 → 报 0.59 cm。
        #
        # 而且"取重叠最多的一组"还会**惩罚延续到松手之后的那段** —— 正常滑动常常
        # 滑过松手时刻，它在窗内的重叠反而小于一段完全落在窗内的短促抖动。
        # 实测 122416：1.19 s 那段（真正的一次滑动，延伸到松手后）重叠只有 0.24 s，
        # 输给了 0.29 s 的抖动段（重叠 0.29 s）。选中的是错的那个。
        sel = [g for g in groups if overlap(g) > 0]
        if not sel:
            return {**res, "reason": "按下期间没有检测到运动"
                                    f"（{len(groups)} 段运动都在按下之前或松手之后）"}
        m0 = min(g[0] for g in sel)
        m1 = max(g[1] for g in sel)
        groups.sort(key=lambda g: g[0])
    elif merge_mode == "button":
        # 运动区间 = 用户按住的那一段，一个样本都不靠静止检测。
        # 代价：按住期间如果停顿很久，窗口会长，零偏误差按 ∝T² 涨。
        # 收益：手抖再大也不会把滑动判成静止（见上面的 slide_wobble）。
        sel = [(down_i, up_i)]
        m0, m1 = down_i, up_i
    else:
        # 「取重叠最大的一组」在**所有组都不与按钮窗重叠**时会退化成「取时间最早的一组」——
        # 那一段往往只是用户按下之前挪了一下手机，不是这次要测的滑动。
        # 注意 `up_i <= down_i` 那个守卫挡不住这里：它只保证按钮窗长度为正，
        # 保证不了选中的运动就是这次滑动。按钮窗完全可以又正又与所有运动组零重叠。
        # 实测：把 100444 的 down_i 换成 up_i（DOWN==UP）它会报 0.06 cm 而不是 15.4 cm。
        groups.sort(key=overlap, reverse=True)
        if overlap(groups[0]) <= 0:
            return {**res, "reason": "按下期间没有检测到运动"
                                    f"（{len(groups)} 段运动都在按下之前或松手之后）"}
        m0, m1 = groups[0]
    res.update(motion=(m0, m1), groups=[tuple(g) for g in groups],
               nMerged=len(sel) if merge_mode == "press" else 1,
               motionOverlap=overlap((sel[0] if merge_mode in ("press", "button")
                                      else groups[0])) * dt)

    # ---- 3. 窗口 = 运动段 ± 余量，两端必须落在静止段里 ----
    m = max(1, int(MARGIN_S / dt))
    w0, w1 = m0 - m, m1 + m
    # 向两侧微调，确保端点是静止样本（余量本身就在静止段里，这里只是兜底）
    while w0 > 0 and not flags[w0]:
        w0 -= 1
    while w1 < n - 1 and not flags[w1]:
        w1 += 1
    w0 = max(0, w0)
    w1 = min(n - 1, w1)
    # 「终点锚点没得用」有两种：① 运动一直延续到搜索区末尾（手机没停稳）；
    # ② 算法要把窗口往右放 MARGIN_S 的余量，但数据在 w1=n-1 处就到头了 ——
    #    这时终点锚点会被裁成只剩几个样本甚至 1 个，白噪声平均不掉。
    #    实测 20261011_100555 尾部静止只有 0.047 s，终点锚点退化成 1 个样本，
    #    它给出的 1.01 cm 完全靠不住，却只由平台宽度（60%）间接暴露。
    end_hit_bound = (m1 >= reg_hi - 1) or (w1 >= n - 1)
    res.update(window=(w0, w1), endHitBound=end_hit_bound)

    if not flags[w0] or not flags[w1]:
        # 两端要分开说。写死一句「松手后没停稳」会让用户朝着错误方向修 ——
        # 实测 20261011_100601 / 100608 是**按下之前**手机一直在动（搜索区里没有一个静止样本），
        # 用户再怎么"松手后保持不动"也没用，得先放稳再按。
        parts = []
        if not flags[w0]:
            parts.append("按下之前手机一直在动")
        if not flags[w1]:
            parts.append("松手后手机没停稳")
        return {**res, "reason": "找不到运动两端的静止段（" + "，".join(parts) + "）"}
    if w1 - w0 < max(5, int(0.05 / dt)):
        return {**res, "reason": "有效窗口太短"}

    # 窗口两端的整段静止（ZUPT 锚点用；在 zupt 内部会被裁到窗口内，
    # 所以实际锚点就是余量那一段静止 —— 窗口**不能**因此变长，否则误差按 T² 涨）
    def run_containing(i):
        a = i
        while a > 0 and flags[a - 1]:
            a -= 1
        b = i + 1
        while b < n and flags[b]:
            b += 1
        return (a, b)
    run0 = run_containing(w0)
    run1 = run_containing(w1)
    if run0[0] == run1[0]:
        return {**res, "ok": True, "distanceM": 0.0, "noMotion": True,
                "imuDistM": 0.0, "axisDistM": 0.0, "leverM": 0.0,
                "straightness": 0.0, "straightRms": 0.0,
                "reason": "窗口两端落在同一段静止里"}
    # STATIC_MIN_S 约束的是**整段静止**，不是裁进窗口之后剩下的那点余量。
    # 后者恒等于 MARGIN_S（实测 0.143 s / 15 个样本，见 exp_anchor.py）——
    # 拿它去卡 STATIC_MIN_S=0.20 会把全部测量干掉，那是约束错了对象。
    #
    # 想让有效锚点变长只有一个办法：把窗口往静止段里外扩。试过了，不行 ——
    # 三条有真值的里两条变差（100459 +5.02%→+9.37%、100529 +3.72%→+6.04%，
    # T 从 2.84 s 涨到 5.04 s），偏差随窗口长度单调增长。
    # 结论：静止段延长窗口**不是**白送样本（零偏会漂，线性模型拟合不住），
    # ∝T² 那一派赢了，√N 那一派输了。所以锚点就这么长，认了。
    #
    # 那 STATIC_MIN_S 还剩下什么用？挡「假静止」：0.2 s 的滚动标准差
    # 根本判不出只有 0.038 s 的静止，而 100555 右侧静止段只有 4 个样本，
    # 它给出的 1.01 cm（平台 0.90~1.51）纯属噪声。整段静止够长，
    # 才说明这里真的静过，值得当锚点。
    if min(run0[1] - run0[0], run1[1] - run1[0]) < min_run:
        s0 = (run0[1] - run0[0]) * dt
        s1 = (run1[1] - run1[0]) * dt
        which = "按下之前" if s0 <= s1 else "松手之后"
        shorter = min(s0, s1)
        return {**res, "reason": f"{which}的静止段只有 {shorter:.2f} s"
                                f"（要 ≥{STATIC_MIN_S:.2f} s 才敢当锚点）—— 停稳后再松手"}
    res.update(run0=run0, run1=run1,
               anchorRunS=((run0[1] - run0[0]) * dt, (run1[1] - run1[0]) * dt))

    # ---- 4. 校准段 = 运动之前那段静止（就是用户「按住后停顿」的那一段）----
    cal = (max(run0[0], m0 - int(CAL_MAX_S / dt)), m0)
    if (cal[1] - cal[0]) * dt < CAL_MIN_S:
        return {**res, "reason": f"运动之前的静止段太短（{(cal[1]-cal[0])*dt:.2f} s）"}
    ca, cb = cal
    res["cal"] = cal
    res["calS"] = (cb - ca) * dt
    gb2 = tuple(statistics.fmean(gyr[i][a] for i in range(ca, cb)) for a in range(3))
    am = tuple(statistics.fmean(acc[i][a] for i in range(ca, cb)) for a in range(3))
    g_mag = vnorm(am)
    # |g| 是个**常量**（当地重力 × 器件比例因子），不该每次从一小段数据里自举：
    # 校准段只有 1~2 s，而手抖是 1~2 Hz 的低频 —— 抖动在这个窗口上根本平均不掉。
    # 实测（仿真 slide_wobble，真值 9.81）：自举值在 9.753~9.794 之间乱跳，
    # 差 0.057 m/s² 的重力残差在 2 s 窗口上积分两步就是 ~11 cm 的假位移。
    # 给了 g_mag_fix（离线标定值）就用它，自举值只留作交叉校验（见 gMagDelta）。
    g_mag_data = g_mag
    if g_mag_fix is not None:
        g_mag = float(g_mag_fix)
    g_ref = vunit(am)
    sigma = max(statistics.pstdev([acc[i][a] for i in range(ca, cb)]) for a in range(3))
    res.update(gMag=g_mag, gMagData=g_mag_data, gMagDelta=g_mag_data - g_mag,
               sigma=sigma, calStaticFrac=sum(
        1 for i in range(ca, cb) if flags[i]) / (cb - ca))
    if not (0.5 * G < g_mag < 1.5 * G):
        return {**res, "reason": f"校准段重力读数异常 |g|={g_mag:.3f}"}
    if res["calStaticFrac"] < CAL_STATIC_FRAC:
        return {**res, "reason": f"校准段不够静止（静止占比 {res['calStaticFrac']*100:.0f}%）"}

    # ---- 5. 姿态（修正参考系）----
    # 按钮窗内一律不做加速度修正：那段时间手机在动，acc 方向不是重力方向。
    # 这条只影响姿态，不影响锚点搜索 —— 锚点仍然用原始终止判据找。
    gate_flags = ([flags[i] and not (down_i <= i <= up_i) for i in range(n)]
                  if gate_motion else flags)
    qs = solve_attitude_v2(acc, gyr, gb2, g_ref, flags, dt, use_flags=gate_flags,
                           i0=ca, g_mag=g_mag)
    a_lin = to_linear(acc, qs, g_mag)

    # ---- 6. 双锚点 ZUPT ----
    dist, vc, p, anchor_drift, run_spread = zupt_two_anchor(a_lin, w0, w1, run0, run1, dt, quad)
    bare = vnorm(integrate_from(a_lin, w0, w1, dt)[1][-1])

    # ---- 6b. 测量头几何（刚体杠杆臂）----
    # 加速度计装在机身某处，用户想测的是机身另一个点。两者**不是同一个点**：
    # 手机一转，即使测量头钉着不动，IMU 也会沿弧线走一段。
    #
    # 刚体关系在**位置层面**是精确等式，不需要角加速度：
    #     p_H(t) = p_I(t) + R(t)·r      ⇒   Δp_H = Δp_I + [R(t1) − R(t0)]·r
    # 这一点很关键：若在加速度层面补偿（a_H = a_I + R(α×r + ω×(ω×r))），
    # 就得对带噪的陀螺读数做数值微分求 α —— 噪声被放大，得不偿失。
    # 位置层面只用到两端姿态，白噪声只影响它们各自的估计误差，不放大。
    #
    # 量级（r=9 cm）：转 30° 时 |ΔR·r| = 2·r·sin(15°) ≈ 4.7 cm。
    # 仿真实测：测量头不动、手机绕它转 30°，v3 报 4.66 cm（真值 0），修正后应回到 0。
    r = tuple(float(x) for x in (head_r or (0.0, 0.0, 0.0)))
    d_lever = vsub(qrot(qs[w1], r), qrot(qs[w0], r))   # [R(t1)−R(t0)]·r，导航系
    p_head = vadd(p[-1], d_lever)
    dist_head = vnorm(p_head)

    # ---- 6c. 沿轨迹主方向的净位移 ----
    # 用户问的是"沿桌边滑了多远"（一个标量），不是"IMU 在三维里飘了多远的模长"。
    # 两者在误差面前不是一回事：模长会把垂直于测量方向的漂移也算进去
    # （z 方向漂 2 cm，15 cm 就报成 15.13 cm），而沿主方向的投影不会。
    # 前提：测量头确实沿一条直线走 —— 用 straightness 把这个前提暴露出来，
    # 轨迹明显不是直线时（比如画弧、抬起来了）投影会低估，诊断会喊。
    p_traj = [vadd(p[k], vsub(qrot(qs[min(w0 + k, n - 1)], r), qrot(qs[w0], r)))
              for k in range(len(p))]
    e = principal_dir(p_traj)
    dist_axis = abs(sum(e[i] * p_head[i] for i in range(3)))
    mu = tuple(sum(q[i] for q in p_traj) / len(p_traj) for i in range(3))
    off = []
    for q_ in p_traj:
        d = vsub(q_, mu)
        t = sum(d[i] * e[i] for i in range(3))
        off.append(vnorm(vsub(d, vscale(e, t))))
    straightness = (max(off) / dist_axis) if dist_axis > 1e-9 else 0.0
    straight_rms = (math.sqrt(sum(x * x for x in off) / len(off)) / dist_axis
                    if dist_axis > 1e-9 else 0.0)

    # 「裸积分」与「去趋势」是残差加速度的两种极端模型：前者当它是 0，后者当它是 A·t+B。
    # 两者相差超过答案本身 ⇒ 答案是两个大数之差，由**模型假设**决定而不是由数据决定，
    # 这时窄平台毫无意义（平台只反映端点抖动，反映不了模型分歧）。
    # 真机实测的分界很干净：有三条有真值（15 cm）的是 3% / 16% / 46%，都没超；
    # 四条没有真值的是 150% / 260% / 590% / 4587%，全超。见 compare_v1_v2.py。
    bare_head = vnorm(vadd(integrate_from(a_lin, w0, w1, dt)[1][-1], d_lever))
    model_gap = abs(bare_head - dist_head)
    model_sensitive = model_gap > max(dist_head, MODEL_GAP_FLOOR_M)

    # ---- 6c. 起手前静止段的角抖动（端点参考稳不稳）----
    # 名字里刻意没有"零偏"：它量的是**手抖**，不是器件漂移 —— 依据见
    # anchor_gyro_tremor 的 docstring（76 s 静止录制证伪了游走假设）。
    # 这里只算、只报，**不**据此拒答、也不乘 T³ 外推：它按场次分层很清楚
    # （桌面 0.0001~0.0004 vs 手持 0.012~0.020），但只有 3 条手持硬真值，
    # 定一条阈值就等于用真值调参。
    tremor, tremor_tau = anchor_gyro_tremor(gyr, run0, dt)

    # ---- 6d. 锚点可信度：松手之后手机到底停了没有 ----
    # 「松手」不等于「停住」：手离开后手机常常还在滑。而**匀速滑动与静止在加速度计上
    # 完全不可分辨**（它测的是比力，匀速时比力就是重力），静止检测必然把这段判成静止。
    # v3 那条教训（匀速段把一次滑动切碎）当时只用来改**运动段**判据，没意识到它同样
    # 毁掉**终点锚点** —— 而终点锚点正是双锚点里"把漂移钉回零"的那一半。
    #
    # 实测 20261011_143734（真值 15 cm，原报 1.44 cm）：
    #   终点锚点段 0.143 s 内积分速度稳定在 −0.31 m/s，运动段速度峰值才 0.314 m/s；
    #   双锚点照假设把它当零，扣掉 15.42 cm **真实位移**，16.19 → 1.59 cm。
    #   同一条的起点锚点 |V0| = 0.007 m/s（真静止）—— 问题只在终点一侧。
    #
    # 这里**只诊断、不拒答**：判据没有干净的分界线。143725 的锚点速度 0.215 m/s
    # 同样接近它自己的运动峰值 0.238 m/s，却给出 +5% 的正确读数（两个错误互相抵消）。
    # 拿 3~6 条真值去拟合一条分界线，等于用真值调参 —— 本项目不做这件事。
    v_win, _ = integrate_from(a_lin, w0, w1, dt)
    v_pk = (max(vnorm(v_win[k]) for k in range(m0 - w0, m1 - w0 + 1))
            if m1 > m0 else 0.0)
    # 下限 0.05 m/s：运动段速度峰值本身就接近零时（122256 只有 0.003 m/s，几乎是没动），
    # 比值是两个噪声相除，没有意义。
    anchor_suspect = v_pk > 0.05 and anchor_drift > 0.5 * v_pk


    # ---- 7. 平台宽度：窗口端点在静止段内挪动，答案变化多少 ----
    grid0, grid1 = [], []
    for k in range(PLATEAU_GRID):
        f = k / (PLATEAU_GRID - 1)
        grid0.append(int(round(run0[0] + f * (m0 - run0[0]))))
        grid1.append(int(round(m1 + f * (run1[1] - 1 - m1))))
    ds = []
    for a0 in grid0:
        for b1 in grid1:
            if b1 - a0 < 5:
                continue
            # 平台扫描的端点在动，杠杆臂项跟着端点姿态走 —— 它必须一起扫，
            # 否则平台只反映 IMU 那部分，掩盖了姿态误差经由杠杆臂的放大。
            p_imu = zupt_two_anchor(a_lin, a0, b1, run_containing(a0),
                                    run_containing(b1), dt, quad)[2][-1]
            dl = vsub(qrot(qs[b1], r), qrot(qs[a0], r))
            ds.append(vnorm(vadd(p_imu, dl)))
    ds.sort()
    plateau = (ds[0] * 100, ds[-1] * 100) if ds else (0.0, 0.0)

    out = {**res, "ok": True, "distanceM": dist_head, "imuDistM": dist,
           "leverM": vnorm(d_lever), "bareM": bare_head, "bareImuM": bare,
           "plateau": plateau,
           "anchorDrift": anchor_drift, "runSpread": run_spread,
           "vPeakMps": v_pk, "anchorSuspect": anchor_suspect,
           "modelGapM": model_gap, "modelSensitive": model_sensitive,
           "T": (w1 - w0) * dt, "holdS": (up_i - down_i) * dt,
           "motionS": (m1 - m0) * dt, "dt": dt,
           # 松手之后到底留了多少静止当终点锚点（秒）。< MARGIN_S 就是"没得用"。
           "tailStaticS": (run1[1] - m1) * dt,
           "flags": flags, "staticRuns": static_runs, "gyroBias": gb2,
           "gRef": g_ref, "anchorRuns": (run0, run1), "motion": (m0, m1),
           # run0/run1 内部是半开区间 [a,b)；Kotlin 侧报闭区间（与 window/motion 一致），
           # 这里显式给出闭区间形式，好让两边的下标能逐字对照。
           "anchorInclusive": ((run0[0], run0[1] - 1), (run1[0], run1[1] - 1)),
           "endHitBound": end_hit_bound,
           # 辨识杠杆臂用：IMU 位移矢量 + 窗口两端姿态。
           # 有了这两个，Δp_H = Δp_I + [R(t1)−R(t0)]·r 可以事后对任意 r 重算，
           # 不必每换一个 r 就把整条链路再跑一遍。
           "pImuVec": tuple(p[-1]), "q0": qs[w0], "q1": qs[w1],
           "dRr": tuple(d_lever),
           "axisDistM": dist_axis, "axis": e,
           "straightness": straightness, "straightRms": straight_rms,
           # 起手前静止段的角抖动（端点参考稳不稳）与它所在的 τ。
           # 只作诊断，不拒答 —— 见 6c 与 anchor_gyro_tremor。
           "startTremorRadS": tremor, "startTremorTauS": tremor_tau}
    warns = []
    if end_hit_bound:
        warns.append("松手后手机在 %0.1f s 内没有停稳（运动一直延续到数据末尾）" % tail_s)
    if out["tailStaticS"] < MARGIN_S:
        warns.append("松手后只留了 %.2f s 静止当终点锚点（引擎要 ≥%.2f s），这个数不可靠 —— "
                     "松手后手别碰、多等一会儿" % (out["tailStaticS"], MARGIN_S))
    # model_sensitive 不再进 warns：它在 §8 直接变成拒答，reason 里已经说了同一件事，
    # 再挂一条"这个数由残差模型决定"的提醒是自相矛盾（一边说别信、一边把数递出去）。
    if anchor_suspect:
        # 不拒答、只提醒：这条判据没有干净分界线（见 6d 的注释）。
        # 但它是 143734 那个 15 cm → 1.44 cm 的病根，值得直说。
        warns.append("终点锚点可疑：标为「静止」的那段里积分速度还有 %.3f m/s，"
                     "与运动段速度峰值 %.3f m/s 同量级 —— 手离开后手机可能还在匀速滑，"
                     "而匀速和静止在加速度计上分不出来。松手后请把手机按住别动"
                     % (anchor_drift, v_pk))
    if warns:
        out["warn"] = "；".join(warns)

    # ---- 8. 模型敏感 → 拒答（不再"给数 + 点名"）----
    # 原则是"宁可拒答，也不要一个像样的假数"。modelSensitive 的物理含义本来就是
    # "这个数由模型假设决定、不是由数据决定" —— 那它不该给数。
    # 原先照给再点名，实测后果：143734 报 1.44 cm 并附带一句"不可当真"，
    # 但 1.44 cm 是个**体面的**数字，用户没有真值可核对，只会照着读。
    # 全库 32 条给数里有 8 条落进这里，其中**只有 143734 有独立真值，而它是错的**；
    # 另一侧已知正确的 5 条（143719/143725/100444/100459/100529）一条都不在名单里。
    #
    # 放在平台扫描**之后**：这样拒答与给数的字段完全一致（含 plateau），
    # 诊断和被拒的那条本身都还能看。代价是一次白扫，可忽略。
    if model_sensitive:
        return {**out, "ok": False,
                "reason": "这个数由残差模型决定、不是由数据决定（裸积分 %.2f cm、"
                          "去趋势 %.2f cm，修正量是答案的 %.0f%%）—— 平台窄只说明端点稳，"
                          "说明不了模型对。请重滑一次：快一点、起手和收手都停稳"
                          % (bare_head * 100, dist_head * 100,
                             model_gap / max(dist_head, 1e-9) * 100)}
    return out


# ------------------------------------------------------------------ 读真机采集

def read_pairs(path):
    ts, rows = [], []
    with open(path, encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            p = line.strip().split(",")
            if len(p) < 4:
                continue
            ts.append(int(p[0]))
            rows.append((float(p[1]), float(p[2]), float(p[3])))
    return ts, rows


def run_dir(d, verbose=True, merge_mode=MERGE_MODE):
    import re
    acc_ts, acc = read_pairs(os.path.join(d, "acc.csv"))
    gyr_ts, gyr = read_pairs(os.path.join(d, "gyr.csv"))
    meta = open(os.path.join(d, "meta.txt"), encoding="utf-8").read()
    m = re.search(r"DOWN\s+event=\d+ ms\s+recv_uptime=\d+ ms\s+recv_elapsed=(\d+) ns", meta)
    if not m:
        return None
    dn = int(m.group(1))
    up = int(re.search(r"UP\s+event=\d+ ms\s+recv_uptime=\d+ ms\s+recv_elapsed=(\d+) ns", meta).group(1))
    n = min(len(acc), len(gyr))
    acc, gyr = acc[:n], gyr[:n]
    dt = statistics.median([acc_ts[i + 1] - acc_ts[i] for i in range(n - 1)]) / 1e9

    def locate(t):
        lo, hi = 0, n - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if acc_ts[mid] < t:
                lo = mid + 1
            else:
                hi = mid
        return lo
    r = measure(acc, gyr, locate(dn), locate(up), dt=dt, merge_mode=merge_mode)
    if verbose:
        nm = os.path.basename(d)
        if not r["ok"]:
            print(f"{nm:>16}  拒答：{r['reason']}")
        else:
            lo, hi = r.get("window", (-1, -1))
            print(f"{nm:>16}  T={r.get('T',0):5.2f}s 按={r.get('holdS',0):5.2f}s "
                  f"动={r.get('motionS',0):5.2f}s "
                  f"窗口#{lo}→#{hi} 静止#{r.get('run0')}{r.get('run1')} "
                  f"|g|={r.get('gMag', 0):.4f} σ_a={r.get('sigma', 0):.4f} "
                  f"★{r['distanceM']*100:7.2f}cm 平台 {r['plateau'][0]:.2f}~{r['plateau'][1]:.2f}"
                  + (f"  〔{r['reason']}〕" if r.get("reason") else "")
                  + ("  ⚠" + r["warn"] if r.get("warn") else ""))
    return r


def main():
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "android", "captures", "slideruler")
    for nm in sorted(os.listdir(root)):
        d = os.path.join(root, nm)
        if os.path.isdir(d):
            run_dir(d)


if __name__ == "__main__":
    main()
