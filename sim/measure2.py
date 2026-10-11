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
                       vcross, vnorm, vscale, vsub, vunit)

# ------------------------------------------------------------------ 参数

STATIC_WIN_S = 0.20      # 静止检测的滑动窗长（边界要锐）
SLOW_WIN_S = 0.80        # 慢动作检测的滑动窗长（白噪声的 σ 不随窗长变，慢信号却会涨）
STATIC_FACTOR = 4.0      # 静止阈值 = 本机噪声地板 × 该系数
STATIC_FLOOR = 0.010     # 噪声地板的硬下限（m/s²），防止极端记录把阈值压到噪声里
STATIC_MIN_S = 0.20      # 一段静止至少持续这么久才算数
GYRO_FLOOR = 0.010       # 陀螺静止阈值的硬下限（rad/s）
CAL_MIN_S = 0.25         # 校准段最短时长（0.25 s 足够：陀螺零偏估计误差 σ/√N ≈ 1.4e-4 rad/s）
CAL_MAX_S = 2.0          # 校准段最长（只用运动前这一段的末尾，别拿太旧的数据）
MARGIN_S = 0.15          # 运动段两侧各留的静止余量（落在静止段内，所以无害）
TAIL_S = 2.5             # UP 之后最多再看这么久，用来找运动结束后的静止段
PRE_S = 2.0              # DOWN 之前最多往回看这么久，用来找校准段
GAP_S = 0.15             # 两段运动之间隔这么久以内的静止，算同一段运动
ALPHA = 0.05             # 互补滤波系数
GATE_RAD = 0.30          # 加速度修正的门限（预测/实测重力方向夹角）
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


def static_flags(acc, gyr, gb, dt):
    """静止标记 + 用到的阈值（全部来源可追溯）。"""
    sd = rolling_max_std(acc, dt)
    floor = percentile(sd, 0.05)
    thr_a = max(STATIC_FACTOR * floor, STATIC_FLOOR)
    om = [vsub(gyr[i], gb) for i in range(len(gyr))]
    on = rolling_max_norm(om, dt)
    ofloor = percentile(on, 0.05)
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

def solve_attitude_v2(acc, gyr, gb, g_ref_body, flags, dt, alpha=ALPHA, gate=GATE_RAD):
    """互补滤波。**加速度修正在静止段才开**，运动中纯陀螺外推。

    v1 的参考是 `g_ref_body`（校准时刻的体重力方向，体坐标里固定不动）——
    等于让滤波器主动抵抗手机的真实旋转：持续角速率 ω 下稳态姿态误差 ≈ ω·dt/α
    （α=0.05、dt=9.5 ms 时 ≈ 0.19·ω），ω=5°/s 就有 1° 误差、重力泄漏 0.17 m/s²。
    v2 换成当前**实测**加速度方向（标准 Mahony），手机转它就跟着转，没有这个稳态误差。
    """
    q = q_from_two_vectors(g_ref_body, (0.0, 0.0, 1.0))
    qs = [q]
    for i in range(len(acc) - 1):
        q_pred = qnorm(qmul(q, q_delta(vsub(gyr[i], gb), dt)))
        if flags[i]:
            e = vcross(qrotT(q_pred, (0.0, 0.0, 1.0)), vunit(acc[i]))
            q = (q_pred if vnorm(e) > gate
                 else qnorm(qmul(q_pred, q_delta(vscale(e, -alpha), dt))))
        else:
            q = q_pred
        qs.append(q)
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


# ------------------------------------------------------------------ 主流程

def measure(acc, gyr, down_i, up_i, dt=None, tail_s=TAIL_S, verbose=False, quad=QUAD):
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
    flags, floor, thr_a, ofloor, thr_g = static_flags(acc, gyr, gb, dt)
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

    # 取与按钮窗重叠最多的那一组
    def overlap(g):
        return max(0, min(g[1], up_i) - max(g[0], down_i))
    groups.sort(key=overlap, reverse=True)
    m0, m1 = groups[0]
    res.update(motion=(m0, m1), groups=[tuple(g) for g in groups])

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
                "reason": "窗口两端落在同一段静止里"}
    res.update(run0=run0, run1=run1)

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
    g_ref = vunit(am)
    sigma = max(statistics.pstdev([acc[i][a] for i in range(ca, cb)]) for a in range(3))
    res.update(gMag=g_mag, sigma=sigma, calStaticFrac=sum(
        1 for i in range(ca, cb) if flags[i]) / (cb - ca))
    if not (0.5 * G < g_mag < 1.5 * G):
        return {**res, "reason": f"校准段重力读数异常 |g|={g_mag:.3f}"}
    if res["calStaticFrac"] < 0.95:
        return {**res, "reason": f"校准段不够静止（静止占比 {res['calStaticFrac']*100:.0f}%）"}

    # ---- 5. 姿态（修正参考系）----
    qs = solve_attitude_v2(acc, gyr, gb2, g_ref, flags, dt)
    a_lin = to_linear(acc, qs, g_mag)

    # ---- 6. 双锚点 ZUPT ----
    dist, vc, p, anchor_drift, run_spread = zupt_two_anchor(a_lin, w0, w1, run0, run1, dt, quad)
    bare = vnorm(integrate_from(a_lin, w0, w1, dt)[1][-1])

    # 「裸积分」与「去趋势」是残差加速度的两种极端模型：前者当它是 0，后者当它是 A·t+B。
    # 两者相差超过答案本身 ⇒ 答案是两个大数之差，由**模型假设**决定而不是由数据决定，
    # 这时窄平台毫无意义（平台只反映端点抖动，反映不了模型分歧）。
    # 真机实测的分界很干净：有三条有真值（15 cm）的是 3% / 16% / 46%，都没超；
    # 四条没有真值的是 150% / 260% / 590% / 4587%，全超。见 compare_v1_v2.py。
    model_gap = abs(bare - dist)
    model_sensitive = model_gap > max(dist, MODEL_GAP_FLOOR_M)

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
            ds.append(zupt_two_anchor(a_lin, a0, b1, run_containing(a0),
                                      run_containing(b1), dt, quad)[0])
    ds.sort()
    plateau = (ds[0] * 100, ds[-1] * 100) if ds else (0.0, 0.0)

    out = {**res, "ok": True, "distanceM": dist, "bareM": bare, "plateau": plateau,
           "anchorDrift": anchor_drift, "runSpread": run_spread,
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
           "endHitBound": end_hit_bound}
    warns = []
    if end_hit_bound:
        warns.append("松手后手机在 %0.1f s 内没有停稳（运动一直延续到数据末尾）" % tail_s)
    if out["tailStaticS"] < MARGIN_S:
        warns.append("松手后只留了 %.2f s 静止当终点锚点（引擎要 ≥%.2f s），这个数不可靠 —— "
                     "松手后手别碰、多等一会儿" % (out["tailStaticS"], MARGIN_S))
    if model_sensitive:
        warns.append("这个数由残差模型决定，不是由数据决定：裸积分 %.2f cm、去趋势 %.2f cm，"
                     "修正量是答案的 %.0f%% —— 平台窄只说明端点稳，说明不了模型对"
                     % (bare * 100, dist * 100, model_gap / max(dist, 1e-9) * 100))
    if warns:
        out["warn"] = "；".join(warns)
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


def run_dir(d, verbose=True):
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
    r = measure(acc, gyr, locate(dn), locate(up), dt=dt)
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
