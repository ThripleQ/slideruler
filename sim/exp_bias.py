"""实验 15：陀螺零偏不是常量 —— 用**所有静态段**估计它在窗口内的游走。

背景（三条 15 cm 硬真值，2026-10-11）：
  · 143734 裸积分 15.78 cm（真值 15）却被去趋势砍到 1.59；
  · 143719 裸积分 4.16 → 去趋势 11.85；
  · 143725 裸积分 20.91 → 去趋势 16.23。
  三条的末端静止段里 a_lin 残余都不为零、且方向一致 ⇒ 姿态漂移是主项。
  而姿态漂移的主因是**陀螺零偏在窗口内游走**（in-run instability ~0.02 rad/s），
  当前实现只用**校准段**一个常值零偏，等于用一个过期 1 秒的值去外推整段。

三种零偏模型对比：
  A 常值：mean(gyr in cal)                        —— 现状
  B 线性插值：两个静止段各自的均值之间线性过渡     —— 用上第二个静止段
  C 最小二乘：所有静止段样本拟合 b(t)=b0+b1·t      —— 用上所有静止段

判据：静止段（未参与估计的那一段！）里 a_lin 的残余模长应当接近 0；
      以及三条硬真值的误差。
"""
import sys
import os
import math
import statistics

sys.path.insert(0, '.')
from exp_real import load, ROOT
from pipeline import (vnorm, vunit, vsub, vscale, vadd, qrotT, vcross,
                      q_from_two_vectors, qnorm, qmul, q_delta)
import measure2 as M

TRUTH = {"20261011_143719": 15.0, "20261011_143725": 15.0, "20261011_143734": 15.0}


def static_runs_of(flags):
    runs = []
    a = None
    for i, f in enumerate(flags):
        if f and a is None:
            a = i
        elif not f and a is not None:
            runs.append((a, i))
            a = None
    if a is not None:
        runs.append((a, len(flags)))
    return [(r[0], r[1]) for r in runs if r[1] - r[0] >= 8]


def bias_constant(gyr, runs, cal):
    ca, cb = cal
    return [statistics.fmean(gyr[i][x] for i in range(ca, cb)) for x in range(3)], None


def bias_linear(gyr, runs, cal):
    """两个静止段（校准段 + 末段）的均值之间线性过渡。返回 (b0, slope 每秒)。"""
    ca, cb = cal
    ends = [r for r in runs if r[1] <= cb + 4]
    if not ends:
        return bias_constant(gyr, runs, cal)
    first = ends[0]
    after = [r for r in runs if r[0] >= cb]
    if not after:
        return bias_constant(gyr, runs, cal)
    last = after[-1]
    t0 = (first[0] + first[1]) / 2.0
    t1 = (last[0] + last[1]) / 2.0
    if t1 - t0 < 10:
        return bias_constant(gyr, runs, cal)
    b0 = [statistics.fmean(gyr[i][x] for i in range(*first)) for x in range(3)]
    b1 = [statistics.fmean(gyr[i][x] for i in range(*last)) for x in range(3)]
    slope = [(b1[x] - b0[x]) / (t1 - t0) for x in range(3)]
    return b0, (slope, t0)


def bias_ls(gyr, runs, cal):
    """所有静止段样本拟合 b(t) = b0 + b1·t（每轴独立最小二乘）。"""
    ts, ys = [], [[] for _ in range(3)]
    for (a, b) in runs:
        for i in range(a, b):
            ts.append(float(i))
            for x in range(3):
                ys[x].append(gyr[i][x])
    n = len(ts)
    if n < 8:
        return bias_constant(gyr, runs, cal)
    mt = sum(ts) / n
    stt = sum((t - mt) ** 2 for t in ts)
    if stt < 1e-9:
        return bias_constant(gyr, runs, cal)
    b0, b1 = [], []
    for x in range(3):
        my = sum(ys[x]) / n
        sty = sum((ts[k] - mt) * (ys[x][k] - my) for k in range(n))
        k1 = sty / stt
        b1.append(k1)
        b0.append(my - k1 * mt)
    return b0, (b1, 0.0)


def attitude(acc, gyr, model, runs, cal, flags, dt, g_mag):
    """用给定零偏模型（可时变）跑互补滤波。返回 qs 与 a_lin。"""
    a_mean = [statistics.fmean(acc[i][x] for i in range(*cal)) for x in range(3)]
    q = q_from_two_vectors(vunit(a_mean), (0.0, 0.0, 1.0))
    qs = [q] * max(0, min(cal[0], len(acc)))
    i0 = max(0, min(cal[0], len(acc) - 1))
    for i in range(max(0, min(i0, len(acc) - 1)), len(acc) - 1):
        b = model(i)
        qp = qnorm(qmul(q, q_delta(vsub(gyr[i], b), dt)))
        # 门控用模长比（F2），静止段才修正
        if flags[i] and abs(vnorm(acc[i]) - g_mag) <= M.GATE_MAG_FRAC * g_mag:
            e = vcross(qrotT(qp, (0.0, 0.0, 1.0)), vunit(acc[i]))
            q = qnorm(qmul(qp, q_delta(vscale(e, -M.ALPHA), dt)))
        else:
            q = qp
        qs.append(q)
    while len(qs) < len(acc):
        qs.append(qs[-1])
    a_lin = M.to_linear(acc, qs, g_mag)
    return qs, a_lin


def residual_in_runs(a_lin, runs):
    out = []
    for (a, b) in runs:
        m = [statistics.fmean(a_lin[i][x] for i in range(a, b)) for x in range(3)]
        out.append((vnorm(m), (a, b)))
    return out


def run_one(nm, verbose=True):
    acc, gyr, dn, up, dt = load(os.path.join(ROOT, nm))
    m = M.measure(acc, gyr, dn, up, dt=dt, merge_mode="press", head_r=None, gate_motion=False)
    if not m["ok"]:
        return None
    flags = m["flags"]
    runs = static_runs_of(flags)
    cal = m["cal"]
    g_mag = m["gMag"]
    w0, w1 = m["window"]
    run0, run1 = m["run0"], m["run1"]

    models = {}
    b0c, _ = bias_constant(gyr, runs, cal)
    models["A 常值"] = (lambda i, b=b0c: b)
    b0l, sl = bias_linear(gyr, runs, cal)
    if sl:
        slope, t0 = sl
        models["B 线性插值"] = (lambda i, b=b0l, s=slope, t=t0:
                           [b[x] + s[x] * (i - t) * dt for x in range(3)])
    else:
        models["B 线性插值"] = models["A 常值"]
    b0s, sl2 = bias_ls(gyr, runs, cal)
    slope2, _ = sl2
    models["C 所有静止段LS"] = (lambda i, b=b0s, s=slope2:
                          [b[x] + s[x] * i * dt for x in range(3)])

    if verbose:
        print("=" * 96)
        print("%s  真值 %.0f cm   静止段 %s   校准段 %s" % (nm, TRUTH[nm], runs, cal))
        print("%-14s %10s %10s %10s %10s %10s" % ("零偏模型", "裸积分cm", "去趋势cm",
                                                  "锚0残余", "锚1残余", "窗口内残余"))
    res = {}
    for tag, model in models.items():
        qs, a_lin = attitude(acc, gyr, model, runs, cal, flags, dt, g_mag)
        v, p = M.integrate_from(a_lin, w0, w1, dt)
        bare = vnorm(p[-1]) * 100
        det = M.zupt_two_anchor(a_lin, w0, w1, run0, run1, dt)[0] * 100
        rr = dict(residual_in_runs(a_lin, runs))
        # 窗口内 a_lin 的均值模长（不参与估计的窗口）
        wm = [statistics.fmean(a_lin[i][x] for i in range(w0, w1 + 1)) for x in range(3)]
        if verbose:
            print("%-14s %10.2f %10.2f %10.4f %10.4f %10.4f" % (
                tag, bare, det, rr.get(runs[0], (0,))[0],
                rr.get(runs[-1], (0,))[0], vnorm(wm)))
        res[tag] = (bare, det)
    return res


def main():
    names = [n for n in TRUTH if os.path.isdir(os.path.join(ROOT, n))]
    for nm in names:
        run_one(nm)
    print("=" * 96)
    print("汇总（三条硬真值的相对误差）")
    print("%-14s %10s %10s %10s" % ("零偏模型", "143734 −89%", "143719 −21%", "143725 +8%"))
    agg = {}
    for nm in names:
        r = run_one(nm, verbose=False)
        for tag, (bare, det) in r.items():
            agg.setdefault(tag, []).append(det)
    for tag, vals in agg.items():
        print("%-14s %s  平均%+.1f%%" % (tag, "  ".join("%9.2f" % v for v in vals),
                                      statistics.fmean((v - 15.0) / 15.0 * 100 for v in vals)))


if __name__ == "__main__":
    main()
