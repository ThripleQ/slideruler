"""窗口端点扫描 + 逐段信号体检（用来判断某次真机采集的距离数字是否可信）

用法: python scan_window.py <采集目录>

它复用 slide_pipeline 的对齐/校准逻辑，然后：
  A. 画出 [起点 x 终点] 距离矩阵 —— 找"平台"。没有平台 = 端点没落在静止段 = 数字不可信。
  B. 逐 0.5 s 打印 |a_lin| 均值/峰值、陀螺范数 —— 看"到底哪段在动"。
  C. 打印按钮窗两端的速度 —— 松手时速度未归零 = 没停稳就松手。
"""

import os
import statistics
import sys

import pipeline as P
from pipeline import G, DeviceProfile, Params, vdot, vnorm, vscale, vunit
import slide_pipeline as S

CAPTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "android", "captures", "slideruler")


def prepare(d):
    acc_ts, acc = S.read_pairs(os.path.join(d, "acc.csv"))
    gyr_ts, gyr = S.read_pairs(os.path.join(d, "gyr.csv"))
    marks = S.read_marks(os.path.join(d, "meta.txt"))
    n = min(len(acc), len(gyr))
    acc, gyr = acc[:n], gyr[:n]
    dts = [acc_ts[i + 1] - acc_ts[i] for i in range(min(len(acc_ts), n) - 1)]
    dt = statistics.median(dts) / 1e9
    t0 = acc_ts[0]
    down_ns = marks["DOWN"]["recv_ns"]
    up_ns = marks["UP"]["recv_ns"]
    lo = S.locate(acc_ts, down_ns)
    hi = S.locate(acc_ts, up_ns)
    cal_n = max(20, min(int(S.CAL_MAX_S / dt), lo))
    cal_lo = lo - cal_n
    b_g = tuple(statistics.fmean(gyr[i][a] for i in range(cal_lo, lo)) for a in range(3))
    a_mean = tuple(statistics.fmean(acc[i][a] for i in range(cal_lo, lo)) for a in range(3))
    a_sig = max(statistics.pstdev([acc[i][a] for i in range(cal_lo, lo)]) for a in range(3))
    g_sig = max(statistics.pstdev([gyr[i][a] for i in range(cal_lo, lo)]) for a in range(3))
    prof = DeviceProfile(dt, b_g, vunit(a_mean), vnorm(a_mean), a_sig, g_sig, cal_n)
    p = Params(static_mode="adaptive", assumed_dt=dt)
    truth = {"N": n - 1, "acc": acc, "gyr": gyr, "dt": dt, "n_pre": lo, "n_mot": hi - lo}
    qs, _ = P.solve_attitude(truth, prof.gyro_bias, prof.gravity_ref, p, dt=dt)
    a_lin = P.remove_gravity(truth, qs, g_mag=G)
    return dict(acc=acc, gyr=gyr, a_lin=a_lin, dt=dt, lo=lo, hi=hi, n=n,
                marks=marks, prof=prof, p=p, truth=truth, t0=t0, acc_ts=acc_ts)


def axis_and_a1(st):
    a_lin, lo, hi = st["a_lin"], st["lo"], st["hi"]
    u, lam1, lam2 = P.pca_axis(a_lin[lo:hi + 1])
    axis = u if vdot(u, (1.0, 0.0, 0.0)) >= 0 else vscale(u, -1.0)
    a1 = [vdot(a_lin[i], axis) for i in range(len(a_lin))]
    return a1, (lam1 / lam2 if lam2 > 0 else float("inf"))


def dist_of(a1, lo, hi, dt):
    v_raw, p_raw = P.integrate(a1, lo, hi, dt=dt)
    T = (hi - lo) * dt
    da = v_raw[-1] / T if T > 1e-9 else 0.0
    a1c = [a1[i] - da if lo <= i <= hi else a1[i] for i in range(len(a1))]
    _, p_z = P.integrate(a1c, lo, hi, dt=dt)
    return abs(p_z[-1]), v_raw[-1]


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else sorted(
        os.path.join(CAPTURES, x) for x in os.listdir(CAPTURES)
        if os.path.isdir(os.path.join(CAPTURES, x)))[-1]
    st = prepare(d)
    dt, lo, hi, n = st["dt"], st["lo"], st["hi"], st["n"]
    a1, lr = axis_and_a1(st)
    print("=" * 96)
    print(f"窗口扫描  {os.path.basename(d)}   dt={dt*1000:.4f}ms  DOWN=#{lo}  UP=#{hi}  λ1/λ2={lr:.1f}")
    print("=" * 96)

    d0, v0 = dist_of(a1, lo, hi, dt)
    print(f"按钮窗 [{lo},{hi}] = {d0*100:.2f} cm    窗口终点原始速度 v_end = {v0:+.4f} m/s "
          f"({'✓ 已停稳' if abs(v0) < 0.05 else '✗ 松手时还在动'})")

    # ---- A. 端点扫描 ----
    starts = sorted(set([max(0, lo - 200), max(0, lo - 100), lo, lo + 60, lo + 120, lo + 240]))
    ends = sorted(set([hi - 240, hi - 120, hi - 60, hi, min(n - 1, hi + 120), min(n - 1, hi + 260)]))
    print()
    print("A · 距离矩阵 (cm)  行=起点样本  列=终点样本   —— 找它有没有『平台』")
    print("      " + "".join(f"{e:>8}" for e in ends))
    for s in starts:
        row = []
        for e in ends:
            if e <= s + 3:
                row.append("   --  ")
            else:
                val, _ = dist_of(a1, s, e, dt)
                row.append(f"{val*100:8.2f}")
        tag = "  ←DOWN" if s == lo else ""
        print(f"  {s:<5}" + "".join(row) + tag)
    print("      " + "".join(f"{'=UP':>8}" if e == hi else "        " for e in ends))
    print("      " + "".join(f"{'↑':>8}" if e == hi else "        " for e in ends))

    # ---- B. 逐 0.5s 信号体检 ----
    print()
    print("B · 逐 0.5 s 信号  (t 相对首样本)")
    w = max(1, int(0.5 / dt))
    print("   时段(s)      |a_lin|均值  |a_lin|峰值   陀螺范数均值   ")
    i = 0
    while i < n:
        j = min(n, i + w)
        am = statistics.fmean(vnorm(st["a_lin"][k]) for k in range(i, j))
        ap = max(vnorm(st["a_lin"][k]) for k in range(i, j))
        gm = statistics.fmean(vnorm(st["gyr"][k]) for k in range(i, j))
        flag = ""
        if i <= lo < j:
            flag = "  ←DOWN"
        if i <= hi < j:
            flag = "  ←UP(松手)"
        bar = "█" * int(min(ap, 2.0) / 0.05)
        print(f"  {i*dt:6.2f}-{j*dt:5.2f}   {am:8.4f}      {ap:8.4f}      {gm:8.5f}   {bar}{flag}")
        i = j

    # ---- C. 检测器对照 ----
    flags = P.detect_static(st["truth"], st["prof"].gyro_bias, st["p"], a_lin=st["a_lin"], profile=st["prof"])
    dlo, dhi = P.motion_window(flags, st["truth"], st["p"])
    print()
    if dlo is None:
        print("C · 静止检测器：没找到运动段")
    else:
        dd, _ = dist_of(a1, dlo, dhi, dt)
        print(f"C · 静止检测器窗 [{dlo},{dhi}] = {dd*100:.2f} cm   （按钮窗 [{lo},{hi}] = {d0*100:.2f} cm）")
    print()


if __name__ == "__main__":
    main()
