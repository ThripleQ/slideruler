"""诊断：双锚点 ZUPT 到底从原始积分里扣掉了什么。

背景（见 docs/v4-head-model.md §5）：
  143734（真值 15 cm）裸积分 15.78 cm ≈ 真值，去趋势后只剩 1.44 cm。
  说明双锚点扣掉的 ~14 cm 几乎全是**真实位移**，不是漂移。
  而 143719 恰好相反：裸积分 4.16 错、去趋势 11.85 对。
  两条数据用同一个模型、判决相反 —— 这是发布算法上最贵的一个洞。

本脚本不做结论，只把中间量摊开：
  · 窗口/运动段/两段锚点裁到窗口内还剩多少
  · 锚点段内**修正前**速度的均值 V0 / V1（这才是"要去掉多少速度"）
  · 拟合出的残差加速度 A 与常值 B
  · 修正项对位移的贡献 ∫c dt（就是被扣掉的那一条）
  · v 与 vc 在窗口内的时间序列（看形状，不看数）

用法：
    python sim/diag_anchor.py                 # 三条硬真值
    python sim/diag_anchor.py 20261011_143734 # 指定一条
"""
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2 as M
from exp_real import load, TRUTH_CM, ROOT


def dump(name, verbose=True):
    r = load(os.path.join(ROOT, name))
    if r is None:
        print("%s 读不了" % name)
        return
    acc, gyr, down_i, up_i, dt = r          # load() 的返回顺序
    m = M.measure(acc, gyr, down_i, up_i, dt=dt, merge_mode="press",
                  head_r=(0.0, -0.09, 0.0))   # 与 exp_real 默认一致
    if not m.get("window"):
        print("%-16s 拒答：%s" % (name, m.get("reason")))
        return

    w0, w1 = m["window"]
    m0, m1 = m["motion"]
    run0, run1 = m["run0"], m["run1"]
    T = (w1 - w0) * dt
    truth = TRUTH_CM.get(name)

    print("=" * 92)
    # 拒答的条也照样诊断（现在 modelSensitive 会直接拒答，而它恰恰最需要看清楚）
    dist_m = m.get("distanceM", m.get("distM"))
    if m.get("ok"):
        verdict = "读数 %.2f cm" % (dist_m * 100)
        if truth:
            verdict += "   误差 %+.1f%%" % ((dist_m * 100 - truth) / truth * 100)
    else:
        verdict = "**拒答**：%s" % (m.get("reason") or "")
    print("%s   真值 %s cm   %s" % (name, truth if truth else "—", verdict))
    print("  窗口 [%d,%d] = %.3f s   运动 [%d,%d] = %.3f s   运动占窗口 %.0f%%"
          % (w0, w1, T, m0, m1, (m1 - m0) * dt, (m1 - m0) / (w1 - w0) * 100))
    print("  完整静止段 run0=%-14s run1=%-14s" % (str(m["run0"]), str(m["run1"])))
    a0 = max(run0[0], w0)
    b0 = min(run0[1], w1 + 1)
    a1 = max(run1[0], w0)
    b1 = min(run1[1], w1 + 1)
    print("  裁到窗口内 run0=[%d,%d] %d 样本 (%.3f s)   run1=[%d,%d] %d 样本 (%.3f s)"
          % (a0, b0, b0 - a0, (b0 - a0) * dt, a1, b1, b1 - a1, (b1 - a1) * dt))

    # --- 复现 measure 内部：取线性加速度、裸速度 ---
    # 注意 measure 的静止标记用的是「最安静 0.5 s」粗估的零偏，**不是**校准段零偏
    # （见 measure2.py §1）。这里必须逐字复现，否则 flags 不同、窗口内的 a_lin 会有
    # 零点几厘米的偏差，诊断出来的 ∫c dt 就跟 measure 报的对不上。
    ca, cb = m["cal"]
    n_all = min(len(acc), len(gyr))
    pre = max(20, min(n_all // 4, int(0.5 / dt)))
    sd0 = M.rolling_max_std(acc, dt)
    i_quiet = min(range(n_all), key=lambda i: sd0[i])
    qa = max(0, i_quiet - pre // 2)
    qb = min(n_all, i_quiet + pre // 2 + 1)
    gb_quiet = tuple(statistics.fmean(gyr[i][a] for i in range(qa, qb)) for a in range(3))
    flags, *_ = M.static_flags(acc, gyr, gb_quiet, dt)
    am = tuple(statistics.fmean(acc[i][a] for i in range(ca, cb)) for a in range(3))
    g_mag = M.vnorm(am)
    gb2 = tuple(statistics.fmean(gyr[i][a] for i in range(ca, cb)) for a in range(3))
    gate_flags = [flags[i] and not (down_i <= i <= up_i) for i in range(n_all)]
    qs = M.solve_attitude_v2(acc, gyr, gb2, M.vunit(am), flags, dt, use_flags=gate_flags,
                             i0=ca, g_mag=g_mag)
    a_lin = M.to_linear(acc, qs, g_mag)
    v, _ = M.integrate_from(a_lin, w0, w1, dt)

    def seg_mean(v_, i0, i1):
        c = i1 - i0
        return [sum(v_[i - w0][x] for i in range(i0, i1)) / c for x in range(3)]

    V0 = seg_mean(v, a0, b0)
    V1 = seg_mean(v, a1, b1)
    t0 = sum((i - w0) * dt for i in range(a0, b0)) / (b0 - a0)
    t1 = sum((i - w0) * dt for i in range(a1, b1)) / (b1 - a1)
    print("  锚点时刻 t0=%.3f s  t1=%.3f s   (t1-t0=%.3f)" % (t0, t1, t1 - t0))
    print("  修正前速度均值 V0 = (%+.4f,%+.4f,%+.4f)  |V0|=%.4f m/s" % (*V0, M.vnorm(V0)))
    print("  修正前速度均值 V1 = (%+.4f,%+.4f,%+.4f)  |V1|=%.4f m/s" % (*V1, M.vnorm(V1)))
    # 滑动自身的速度量级做对照：运动段中点附近的速度峰
    vpk = max(M.vnorm(v[i - w0]) for i in range(m0, m1))
    print("  运动段内速度峰值 %.4f m/s   —— 锚点速度若与它同量级，说明锚点里其实还在动"
          % vpk)

    AB = []
    for x in range(3):
        det = t0 * 1.0 - 1.0 * t1
        AB.append(((V0[x] * 1.0 - 1.0 * V1[x]) / det,
                   (t0 * V1[x] - V0[x] * t1) / det))
    print("  拟合残差加速度 A = (%+.4f,%+.4f,%+.4f)  |A|=%.4f m/s²"
          % (*[AB[x][0] for x in range(3)], M.vnorm([AB[x][0] for x in range(3)])))
    print("  拟合常值      B = (%+.4f,%+.4f,%+.4f)  |B|=%.4f m/s"
          % (*[AB[x][1] for x in range(3)], M.vnorm([AB[x][1] for x in range(3)])))

    # 修正项对位移的贡献 ∫c dt（在窗口 [0,T] 上）
    c_int = [AB[x][0] * T ** 2 / 2 + AB[x][1] * T for x in range(3)]
    print("  被扣掉的位移 ∫c dt = (%+.4f,%+.4f,%+.4f)  |·| = %.2f cm   ← 这就是那 14 cm"
          % (*c_int, M.vnorm(c_int) * 100))

    p_raw = M.integrate_from(a_lin, w0, w1, dt)[1]
    print("  裸积分末端位移 |p_raw| = %.2f cm" % (M.vnorm(p_raw[-1]) * 100))
    print("  去趋势末端位移 |p|     = %.2f cm   （= 裸积分 − 被扣掉的）"
          % (M.vnorm(M.zupt_two_anchor(a_lin, w0, w1, run0, run1, dt)[2][-1]) * 100))

    if verbose:
        print("  窗口内速度序列（每 %d 样本一行，单位 m/s）：" % max(1, (w1 - w0) // 24))
        step = max(1, (w1 - w0) // 24)
        print("    %7s %9s %9s %9s %9s %9s  %s" %
              ("t s", "vx", "vy", "vz", "|v|", "|vc|", "标记"))
        dist, vc, p, ad, sp = M.zupt_two_anchor(a_lin, w0, w1, run0, run1, dt)
        for k in range(0, len(v), step):
            i = w0 + k
            tag = ""
            if m0 <= i <= m1:
                tag += "运动 "
            if a0 <= i < b0 or a1 <= i < b1:
                tag += "锚点"
            print("    %7.3f %+9.4f %+9.4f %+9.4f %9.4f %9.4f  %s"
                  % (k * dt, v[k][0], v[k][1], v[k][2], M.vnorm(v[k]),
                     M.vnorm(vc[k]), tag))


def main():
    names = [a for a in sys.argv[1:] if a.startswith("2026")]
    if not names:
        names = sorted(TRUTH_CM)
    for nm in names:
        dump(nm)


if __name__ == "__main__":
    main()
