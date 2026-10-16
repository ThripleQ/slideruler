"""v5 批量估计（**已否定，未启用**）：把误差修在**姿态**层面，不在位置层面拟合多项式。

> ## 结论：这条路线被仿真器否掉了，**没有同步到 Kotlin**
>
> 想法本身是对的（下面"为什么"那段仍然成立），但**做法**输给了现有双锚点：
>
> | 场景（K40 量级噪声 + 零偏游走 0.02） | v4 双锚点 | v5 批量 |
> |---|---|---|
> | slide T=0.7 s | +3.3%（σ 1.2） | −4.2%（σ 4.3） |
> | slide T=1.0 s | −3.6%（σ 2.8） | −1.7%（σ 9.6） |
> | slide T=1.4 s | **+8.0%**（σ 7.0） | **+47.4%**（σ 13.4） |
> | slide T=2.0 s | +102.7% | +113.8% |
> | static_rotate30（真值 **0**） | 5.9 cm | **13.0 cm** |
>
> σ 大 2~5 倍、长窗口更差、在"真值为 0 的纯旋转"上把虚构量翻倍。
>
> ### 为什么必然输：δp/δθ = ½·g·T²
> 姿态误差到位移的放大系数是 `½·g·T²`：T=1.4 s 时 ≈ **9.7 m/rad**。
> 而静止段的"倾斜观测"噪声实测约 **1°**（受加速度计低频噪声限制，
> 见 `exp_bias.py` 的拟合残差 0.58~1.01°，且不随样本数下降）。
> 拿 1° 的观测去旋转整段姿态 → 直接灌进 5~10 cm。
> 双锚点在**位置空间**扣漂移，对姿态误差的敏感度低得多（线性残差模型天然吸收斜坡），
> 所以它反而稳。**"更物理"不等于"更准" —— 精度由观测噪声决定，不是由模型优雅程度决定。**
>
> ### 顺带否掉的第二个想法：用**所有**静止样本约束 v=0（见 exp_zupt.py ③④）
> 静止标记的样本**在运动中间也会出现**（匀速段与静止在加速度计上不可分辨 —— v3 就发现过，
> 真机 122416 一次按住被切成 5 段）。这些样本的速度**不是零**，是滑动的速度。
> 拿它们去约束"v=0"是错的：仿真里 T=1.4 s 从 +23% 掉到 −52%，T=2.0 s 掉到 −69%。
> **只有"运动开始之前"和"运动结束之后"的静止段才有 v=0。**

## 为什么（想法本身仍然是对的）

v2~v4 的去趋势是「在位置/加速度层面拟合一个自由多项式」。它有两个结构性毛病：

1. **模型不对**。真机的残差加速度不是常量、也不是线性函数，它是
   `a_res(t) = θ(t) × g`，而 θ(t) = 姿态误差 ≈ **陀螺零偏误差 × 时间 + 初值**，
   也就是**斜坡**。用常量+线性去拟合一条斜坡，参数被错配到"常量"上，
   拟合出的修正量与实际残差形状不一致 —— 143734 就是这么被砍掉 14 cm 的
   （裸积分 15.78 cm ≈ 真值 15，去趋势后 1.59 cm）。

2. **只用两个点**。双锚点只用两侧静止段的**均值**，两个方程解两个未知数，
   余量全浪费。而静止段内每一**个**样本都是一条"此刻重力方向应当是竖直"的观测，
   143734 有 190 条静态样本 → 380 个约束，足够把 4 个未知数钉死。

## 怎么做

未知量只有 4 个（倾斜误差在水平面内的 2+2）：
    θ(t) = θ_A + δb·(t − t_A)，   θ_A、δb 都是**水平**矢量（沿竖直的分量不可观测）
观测：每条静止样本 i 给 2 条方程。取 u_i = qrot(q_i, acc_i)/|acc_i|（导航系比力方向，
静止时应当正好是竖直 ẑ），要求修正后 u_i 与 ẑ 重合：
    θ(t_i) × u_i = ẑ − u_i   ⇒   θ_x = −d_y/u_z,  θ_y = +d_x/u_z   （d = ẑ − u）
再对所有静态样本做**线性最小二乘**（θ_x、θ_y 各自对时间回归）→ θ_A、δb。
然后 `q ← exp(θ_A)·q`、`bias ← bias + δb`，重新积分，迭代 2~3 次。
加速度零偏 b_a 由静止段残差的均值直接读出（不再是自由拟合的常数项）。

## 与旧方法的关系

旧方法的「常数项」≈ b_a（可为，但被噪声和错配污染）；
旧方法的「线性项」≈ 把姿态漂移的**位置后果**（∝t³）硬塞进一个 ∝t 的项里 ——
这正是误差最大的地方。新方法让漂移在**源头**（姿态）被扣掉。
"""
import math
import statistics

from pipeline import (qrot, qrotT, qmul, qnorm, q_delta, q_from_two_vectors,
                      vunit, vnorm, vsub, vadd, vscale, vcross)
import measure2 as M

ITER = 3          # 迭代次数（线性化是小角度假设，迭代几次收敛）


def _tilt_obs(acc, qs, i, g_mag=None):
    """第 i 个静止样本给出的**导航系**倾斜误差 (θ_x, θ_y)，以及该样本的姿态 q。

    约定：qrot(q, v) 把**体系**矢量转到导航系；静态时 qrot(q, acc) 应当 = (0,0,|g|)。
    倾斜误差 δ（导航系左乘）满足 q_true = δ ⊗ q̂ ⇒ u_true = u + δ × u。
    要 u_true = ẑ：δ × u = ẑ − u = d。
      δ_y·u_z = d_x,  −δ_x·u_z = d_y  ⇒  δ_x = −d_y/u_z,  δ_y = +d_x/u_z
    沿竖直的 δ_z 不可观测（它不产生重力泄漏）→ 置 0，这是一个规范选择。
    """
    n = qrot(qs[i], acc[i])
    mag = vnorm(n)
    if mag < 1e-9:
        return None
    # 运动中的 acc = 重力 + 运动加速度，方向不是重力方向 —— 用模长比挡掉
    if g_mag and abs(mag - g_mag) > M.GATE_MAG_FRAC * g_mag:
        return None
    u = vscale(n, 1.0 / mag)
    if u[2] < 0.8:            # 姿态错到 37° 以上：线性化失效，丢掉这条
        return None
    d = (0.0 - u[0], 0.0 - u[1], 1.0 - u[2])     # ẑ − u
    return (-d[1] / u[2], d[0] / u[2], 0.0)


def _tilt_obs_body(acc, qs, i, g_mag=None):
    """同上，但把倾斜误差转到**体系**：θ_body = R(q)⁻¹·θ_nav。

    为什么必须转到体系：陀螺零偏误差在体系里是常量，手机一转，导航系里的
    倾斜误差就不再是"一条直线"而是 R(t)·(Δb·t) —— 一条绕着重力轴转的曲线。
    实测证据：`static_rotate30`（手机绕测量头转 30°）上，导航系直线模型报出
    12.8 cm 的纯虚构（真值 0），而体系模型没有这个问题。
    """
    ob = _tilt_obs(acc, qs, i, g_mag)
    if ob is None:
        return None
    return qrotT(qs[i], ob)


def _fit_tilt(acc, qs, static_runs, dt, i_ref, i_min=0, g_mag=None, frame="body"):
    """对所有静止样本做线性回归 θ(t) = A + B·(t−t_ref)。

    frame="nav"  ：直接在导航系里拟合（旧版）—— 只在手机不转时才对
    frame="body" ：先转到体系再拟合 —— 零偏误差在体系里才真正是直线的
    """
    ts, obs = [], []
    for (a, b) in static_runs:
        for i in range(max(a, i_min), b):
            ob = _tilt_obs_body(acc, qs, i, g_mag) if frame == "body" else _tilt_obs(acc, qs, i, g_mag)
            if ob is None:
                continue
            ts.append((i - i_ref) * dt)
            obs.append(ob)
    n = len(ts)
    if n < 12:
        return None
    mt = sum(ts) / n
    stt = sum((t - mt) ** 2 for t in ts)
    out = []
    for c in range(3):
        y = [o[c] for o in obs]
        my = sum(y) / n
        if stt < 1e-12:
            out.append((my, 0.0))
            continue
        sty = sum((ts[k] - mt) * (y[k] - my) for k in range(n))
        b = sty / stt
        out.append((my - b * mt, b))
    # 残差（拟合后仍解释不掉的倾斜）—— 也就是"不可约"的那部分
    res = []
    for k in range(n):
        e = [obs[k][c] - (out[c][0] + out[c][1] * ts[k]) for c in range(3)]
        res.append(vnorm(e))
    return {"a": tuple(o[0] for o in out), "b": tuple(o[1] for o in out),
            "residRad": statistics.fmean(res), "n": n, "frame": frame}


def fit_attitude(acc, gyr, gb, q_init, flags, static_runs, dt, i0, g_mag,
                 iters=ITER, verbose=False, frame="body"):
    """迭代求解姿态倾斜误差 + 陀螺零偏误差。

    返回 (qs, bias, info)。qs[i] 对应 acc[i]，与 measure2.solve_attitude_v2 同约定。

    ⚠️ solve_attitude_v2 只接受"初始重力方向（**体系**）"而不是四元数，
    所以要把它转回去：q0 把 nav-up 映到体系的矢量是 qrotT(q0, ẑ)。
    早先写成 qrot(q0, ẑ)（转反了），会让修正越修越歪 —— 22° 的假倾斜就是这么来的。
    """
    bias = list(gb)
    q0 = q_init
    info = {}
    for k in range(iters):
        body_up = qrotT(q0, (0.0, 0.0, 1.0))   # nav 上方向在体系里的表示
        qs = M.solve_attitude_v2(acc, gyr, tuple(bias), body_up, flags, dt,
                                 i0=i0, g_mag=g_mag)
        fit = _fit_tilt(acc, qs, static_runs, dt, i_ref=i0, i_min=i0,
                        g_mag=g_mag, frame=frame)
        if fit is None:
            info["fit"] = None
            return qs, tuple(bias), info
        a, b = fit["a"], fit["b"]
        info = {"fit": fit, "iter": k + 1, "tiltA": vnorm(a), "tiltSlope": vnorm(b),
                "biasAdd": b}
        if verbose:
            print("    迭代%d: 倾斜初值 %.3f°  斜率 %.4f rad/s  残差 %.3f°"
                  % (k + 1, math.degrees(vnorm(a)), vnorm(b),
                     math.degrees(fit["residRad"])))
        # 更新初值姿态：误差在**体系**里是 a（frame="body"），转到导航系再加到初值姿态上。
        # frame="nav" 时 a 本身就在导航系。
        nav_a = qrot(q0, a) if frame == "body" else a
        q0 = qnorm(qmul(q_delta(nav_a, 1.0), q0))
        # 零偏误差在体系里（frame="body"）才直接就是体系零偏的修正量；
        # nav 模式下只能近似当作体系的 x/y（小倾角时成立）。
        # 符号见 fit_attitude 的推导：观测斜率 = −R·Δb ⇒ Δb = −斜率。
        for c in range(3):
            bias[c] -= b[c]
        if vnorm(a) < 1e-6 and vnorm(b) < 1e-7:
            break
    body_up = qrotT(q0, (0.0, 0.0, 1.0))
    qs = M.solve_attitude_v2(acc, gyr, tuple(bias), body_up, flags, dt,
                             i0=i0, g_mag=g_mag)
    return qs, tuple(bias), info


def measure_batch(acc, gyr, down_i, up_i, dt=None, head_r=None,
                  merge_mode=M.MERGE_MODE, verbose=False, frame="body"):
    """与 measure2.measure 同接口的批量估计版本。

    分段/静止检测/拒答判据**完全复用 measure2.measure**（那部分已被反复验证），
    只替换估计内核：姿态批量对齐 + 加速度零偏直读 + 直接积分（不做多项式去趋势）。
    """
    base = M.measure(acc, gyr, down_i, up_i, dt=dt, merge_mode=merge_mode,
                     head_r=head_r, gate_motion=False)
    if not base.get("ok"):
        return base
    if base.get("noMotion"):
        return {**base, "batch": True}
    dt = base["dt"]
    n = len(acc)
    flags = base["flags"]
    static_runs = [(r[0], r[1]) for r in base["runs"]] if "runs" in base else None
    if not static_runs:
        # measure2 只在返回值里给了 run0/run1；其余段自己再扫一遍
        static_runs = []
        a = None
        for i in range(n):
            if flags[i] and a is None:
                a = i
            elif not flags[i] and a is not None:
                static_runs.append((a, i))
                a = None
        if a is not None:
            static_runs.append((a, n))
        static_runs = [r for r in static_runs if r[1] - r[0] >= 8]

    i0 = base["cal"][0]
    g_mag = base["gMag"]
    a_mean = [statistics.fmean(acc[i][x] for i in range(*base["cal"])) for x in range(3)]
    q_init = q_from_two_vectors(vunit(a_mean), (0.0, 0.0, 1.0))
    gb = [statistics.fmean(gyr[i][x] for i in range(*base["cal"])) for x in range(3)]

    qs, bias, info = fit_attitude(acc, gyr, gb, q_init, flags, static_runs, dt,
                                 i0, g_mag, verbose=verbose, frame=frame)
    a_lin = M.to_linear(acc, qs, g_mag)
    # 加速度零偏：静止段残差的均值（物理直读，不是自由拟合常数项）
    stat_idx = [i for (x, y) in static_runs for i in range(x, y)]
    b_a = [statistics.fmean(a_lin[i][k] for i in stat_idx) for k in range(3)]
    a_corr = [[a_lin[i][k] - b_a[k] for k in range(3)] for i in range(n)]

    w0, w1 = base["window"]
    m0, m1 = base["motion"]
    run0, run1 = base["run0"], base["run1"]
    r = tuple(float(x) for x in (head_r or (0.0, 0.0, 0.0)))
    d_lever = vsub(qrot(qs[w1], r), qrot(qs[w0], r))

    # 从窗口起点积起（那里在静止段内，v≈0），读到窗口止
    v, p = M.integrate_from(a_corr, w0, w1, dt)
    p_full = vadd(p[-1], d_lever)
    # 仍然报"裸积分"以便与旧版对照：同一姿态，但不扣 b_a、不去趋势
    v_bare, p_bare = M.integrate_from(a_lin, w0, w1, dt)
    bare_full = vnorm(vadd(p_bare[-1], d_lever))

    # 静止段内修正后的残余速度（应当≈0）——这是新方法的质量指标
    vmean = {}
    for (x, y) in static_runs:
        idx = [i - w0 for i in range(max(x, w0), min(y, w1 + 1))]
        if len(idx) >= 3:
            vmean[(x, y)] = vnorm([statistics.fmean(v[i][k] for i in idx)
                                   for k in range(3)])

    out = {**base,
           "batch": True,
           "distanceM": (vnorm(p_full) if head_r else vnorm(p[-1])),
           "imuDistM": vnorm(p[-1]),
           "leverM": vnorm(d_lever),
           "bareM": bare_full,
           "modelGapM": abs(bare_full - vnorm(p_full if head_r else p[-1])),
           "staticResidMS": max(vmean.values()) if vmean else 0.0,
           "accBias": tuple(b_a),
           "biasEst": tuple(bias),
           "tiltResidRad": (info.get("fit") or {}).get("residRad", 0.0),
           "tiltA": info.get("tiltA", 0.0),
           "tiltSlope": info.get("tiltSlope", 0.0),
           "qs": qs, "aCorr": a_corr}
    out["modelSensitive"] = out["modelGapM"] > max(out["distanceM"], M.MODEL_GAP_FLOOR_M)
    return out


# ---------------------------------------------------------------- 真机回归

def run_dir(d, **kw):
    import os
    import re
    acc_ts, acc = M.read_pairs(os.path.join(d, "acc.csv"))
    gyr_ts, gyr = M.read_pairs(os.path.join(d, "gyr.csv"))
    meta = open(os.path.join(d, "meta.txt"), encoding="utf-8").read()
    m = re.search(r"DOWN\s+event=\d+ ms\s+recv_uptime=\d+ ms\s+recv_elapsed=(\d+) ns", meta)
    if not m:
        return None
    dn = int(m.group(1))
    up = int(re.search(r"UP\s+event=\d+ ms\s+recv_uptime=\d+ ms\s+recv_elapsed=(\d+) ns",
                       meta).group(1))
    down_i = M.nearest(acc_ts, dn)
    up_i = M.nearest(acc_ts, up)
    return measure_batch(acc, gyr, down_i, up_i, **kw)


if __name__ == "__main__":
    import os
    import sys
    root = os.path.join(os.path.dirname(__file__), "..", "android", "captures", "slideruler")
    truth = {"20261011_143719": 15.0, "20261011_143725": 15.0, "20261011_143734": 15.0}
    print("%-16s %8s %10s %10s %10s %9s" % ("采集", "真值", "v4去趋势", "v5批量", "静止残余", "倾斜残差"))
    for nm in sorted(os.listdir(root)):
        d = os.path.join(root, nm)
        if not os.path.isdir(d):
            continue
        v4 = M.run_dir(d, verbose=False)
        v5 = run_dir(d)
        if v4 is None or v5 is None:
            continue
        t = truth.get(nm)
        def fmt(r):
            if r is None or not r.get("ok"):
                return "拒" if r else "-"
            return "%.2f" % (r["distanceM"] * 100)
        if t or nm in truth:
            print("%-16s %8s %10s %10s %9.3f %8.3f°" % (
                nm, "%.0f" % t if t else "-", fmt(v4), fmt(v5),
                v5.get("staticResidMS", 0.0), math.degrees(v5.get("tiltResidRad", 0.0))))
