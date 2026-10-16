"""实验 17：位置空间的漂移模型该用哪一种？（对照四种 ZUPT 变体）

背景：实验 16 证明"修在姿态层面"输给"修在位置层面"，原因是放大系数
       δp/δθ = ½·g·T²（T=1.4 s 时 ≈9.7 m/rad，1° 观测噪声 → 10 cm）。
而现在的双锚点（位置层面）用的是**线性**残差模型，只用两侧静止段的**均值**：
      v_corr(t) = A·t + B      （a_corr = A，常量）
它对应"残差加速度是常量"。但实测残差加速度是**斜坡**（姿态误差 ∝ t，见 exp_bias.py
里零偏游走 0.006 rad/s 的斜率），斜坡模型应当是
      v_corr(t) = A·t + B·t²/2  （a_corr = A + B·t，线性）
也就是代码里一直没开（且被判定"更差"）的 quad。它当时更差的原因很可能是
**只用两个均值点、6 约束解 6 未知数、毫无余量**（注释里自己写了"矩阵接近奇异"）。
本实验把"用哪些样本"和"用几阶模型"分开，做 2×2 对照：

  ① 均值 + 线性（现状）        ② 均值 + 二次
  ③ 全静止样本 LS + 线性       ④ 全静止样本 LS + 二次

判据：三条 15 cm 硬真值 + 仿真器（真值精确、零偏游走可控）。
"""
import sys
import os
import math
import statistics

sys.path.insert(0, '.')

from pipeline import vnorm, vsub, vscale, vadd
import measure2 as M


def zupt_variant(a_lin, w0, w1, run0, run1, dt, mode="mean", quad=False,
                 static_runs=None):
    """mode="mean"：只用两侧静止段的均值（现状）
       mode="all" ：用窗口内**所有**静止样本做最小二乘
    """
    v, _ = M.integrate_from(a_lin, w0, w1, dt)
    r0 = (run0[0], run0[1])
    r1 = (run1[0], run1[1])

    def idx_of(run):
        return [i for i in range(max(run[0], w0), min(run[1], w1 + 1))]

    if mode == "all" and static_runs:
        idx = [i for i in range(w0, w1 + 1)
               if any(x <= i < y for (x, y) in static_runs)]
        if len(idx) < 6:
            idx = idx_of(r0) + idx_of(r1)
    else:
        idx = idx_of(r0) + idx_of(r1)
    if len(idx) < 4:
        return M.zupt_two_anchor(a_lin, w0, w1, run0, run1, dt, quad)[0], None

    # 设计矩阵：v_corr = A·t + B·t²/2（quad）或 A·t + B（linear）
    rows = [(k * dt, (k * dt) ** 2 / 2.0) if quad else (k * dt, 1.0)
            for k in (i - w0 for i in idx)]
    # 每轴独立 2×2 最小二乘
    s11 = sum(r[0] * r[0] for r in rows)
    s12 = sum(r[0] * r[1] for r in rows)
    s22 = sum(r[1] * r[1] for r in rows)
    det = s11 * s22 - s12 * s12
    AB = []
    for x in range(3):
        y = [v[i - w0][x] for i in idx]
        b1 = sum(rows[j][0] * y[j] for j in range(len(idx)))
        b2 = sum(rows[j][1] * y[j] for j in range(len(idx)))
        if abs(det) < 1e-15:
            AB.append((0.0, 0.0))
        else:
            AB.append(((b1 * s22 - s12 * b2) / det, (s11 * b2 - b1 * s12) / det))

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
    # 静止样本上的残余速度（越小说明模型越贴合"静止段速度该为 0"这个事实）
    resid = statistics.fmean([vnorm(vc[i - w0]) for i in idx])
    return vnorm(p[-1]), resid


VARIANTS = [
    ("① 均值+线性(现状)", dict(mode="mean", quad=False)),
    ("② 均值+二次", dict(mode="mean", quad=True)),
    ("③ 全静止+线性", dict(mode="all", quad=False)),
    ("④ 全静止+二次", dict(mode="all", quad=True)),
]


def static_runs_of(flags, min_len=8):
    runs, a = [], None
    for i, f in enumerate(flags):
        if f and a is None:
            a = i
        elif not f and a is not None:
            runs.append((a, i))
            a = None
    if a is not None:
        runs.append((a, len(flags)))
    return [r for r in runs if r[1] - r[0] >= min_len]


def real_table(names, truth_map):
    from exp_real import load, ROOT
    print("### 真机三条硬真值（15 cm 直尺）")
    print("%-18s %6s | %10s %10s %10s %10s" % ("采集", "真值", *[v[0][2:] for v in VARIANTS]))
    errs = {t: [] for t, _ in VARIANTS}
    for nm in names:
        acc, gyr, dn, up, dt = load(os.path.join(ROOT, nm))
        m = M.measure(acc, gyr, dn, up, dt=dt, merge_mode="press", head_r=None,
                      gate_motion=False)
        if not m["ok"]:
            print("%-18s 拒答" % nm)
            continue
        w0, w1 = m["window"]
        gb2 = tuple(statistics.fmean(gyr[i][a] for i in range(*m["cal"])) for a in range(3))
        am = tuple(statistics.fmean(acc[i][a] for i in range(*m["cal"])) for a in range(3))
        from pipeline import vunit
        qs = M.solve_attitude_v2(acc, gyr, gb2, vunit(am), m["flags"], dt, i0=m["cal"][0],
                                 g_mag=m["gMag"])
        a_lin = M.to_linear(acc, qs, m["gMag"])
        runs = static_runs_of(m["flags"])
        t = truth_map[nm]
        row = []
        for tag, kw in VARIANTS:
            d, _ = zupt_variant(a_lin, w0, w1, m["run0"], m["run1"], dt,
                                static_runs=runs, **kw)
            row.append(d * 100)
            errs[tag].append((d * 100 - t) / t * 100)
        print("%-18s %6.0f | %s   （真值 %.0f）" % (
            nm, t, " ".join("%10.2f" % x for x in row), t))
    print("%-18s %6s | %s" % ("平均相对误差", "", " ".join(
        "%9.1f%%" % statistics.fmean(errs[tag]) for tag, _ in VARIANTS)))
    print("%-18s %6s | %s" % ("最大 |误差|", "", " ".join(
        "%9.1f%%" % max(abs(e) for e in errs[tag]) for tag, _ in VARIANTS)))


def sim_table(cfg_kw, n=8):
    from simgen import SimCfg, scenario, gen, down_i_of
    from pipeline import vunit
    print("### 仿真（真值精确，零偏游走 0.02 rad/s/√s）")
    cases = [("slide_T0.7", "slide", dict(T=5.0, d=0.15, move_s=0.7)),
             ("slide_T1.0", "slide", dict(T=5.0, d=0.15, move_s=1.0)),
             ("slide_T1.4", "slide", dict(T=5.0, d=0.15, move_s=1.4)),
             ("slide_T2.0", "slide", dict(T=5.0, d=0.15, move_s=2.0)),
             ("slide_rot_T1.4", "slide_rotate", dict(T=5.0, d=0.15, move_s=1.4, rot_deg=15))]
    print("%-18s %7s | %s" % ("场景", "真值", " ".join("%12s" % v[0][2:] for v in VARIANTS)))
    for name, sname, kw in cases:
        acc_by = {t: [] for t, _ in VARIANTS}
        for s in range(n):
            spec = scenario(sname, cfg=SimCfg(**cfg_kw), r_body=(0.0, -0.09, 0.0), **kw)
            spec["cfg"].seed = s + 1
            r = gen(spec["cfg"], spec["head_pos"], spec["att"], spec["r_body"], spec["t_end"])
            dt = 1 / spec["cfg"].hz
            dn = down_i_of(r["acc_ts"], spec["tA"] - 0.25)
            up = down_i_of(r["acc_ts"], spec["tB"] + 0.25)
            m = M.measure(r["acc"], r["gyr"], dn, up, dt=dt, merge_mode="press",
                          head_r=(0.0, -0.09, 0.0), gate_motion=False)
            if not m["ok"]:
                continue
            w0, w1 = m["window"]
            gb2 = tuple(statistics.fmean(r["gyr"][i][a] for i in range(*m["cal"]))
                        for a in range(3))
            am = tuple(statistics.fmean(r["acc"][i][a] for i in range(*m["cal"]))
                       for a in range(3))
            qs = M.solve_attitude_v2(r["acc"], r["gyr"], gb2, vunit(am), m["flags"], dt,
                                     i0=m["cal"][0], g_mag=m["gMag"])
            a_lin = M.to_linear(r["acc"], qs, m["gMag"])
            runs = static_runs_of(m["flags"])
            for tag, kwv in VARIANTS:
                d, _ = zupt_variant(a_lin, w0, w1, m["run0"], m["run1"], dt,
                                    static_runs=runs, **kwv)
                acc_by[tag].append(d * 100)
        t = 15.0
        cells = []
        for tag, _ in VARIANTS:
            xs = acc_by[tag]
            if xs:
                cells.append("%5.2f %+5.1f%%" % (statistics.fmean(xs),
                                                 statistics.fmean(xs) - t) if False
                             else "%5.2f%+6.0f%%" % (
                                 statistics.fmean(xs),
                                 (statistics.fmean(xs) - t) / t * 100))
            else:
                cells.append("%12s" % "-")
        print("%-18s %7.0f | %s" % (name, t, " ".join("%12s" % c for c in cells)))


if __name__ == "__main__":
    NOISY = dict(acc_noise=0.010, gyr_noise=0.002, acc_bias=(0.020, -0.010, 0.030),
                 gyr_bias=(0.002, -0.001, 0.003), acc_scale=(0.011, 0.011, 0.011),
                 acc_misalign=0.004, gyr_bias_rw=0.02)
    TRUTH = {"20261011_143719": 15.0, "20261011_143725": 15.0, "20261011_143734": 15.0}
    real_table(sorted(TRUTH), TRUTH)
    print()
    sim_table(NOISY)
