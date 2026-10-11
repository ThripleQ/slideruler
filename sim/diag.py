"""逐采集体检：这条记录到底有没有动、校准段干不干净、距离是怎么来的。

起因：用户报「按着不动也测出很长一段距离」。

用法:
    python diag.py <采集目录> [<采集目录> ...]
    python diag.py --all          # captures/slideruler 下全部
"""
import math
import os
import re
import statistics
import sys

import pipeline as P
from pipeline import G, vdot, vnorm, vunit

HERE = os.path.dirname(os.path.abspath(__file__))
CAP_ROOT = os.path.join(HERE, "..", "android", "captures", "slideruler")


# ----------------------------------------------------------------- 读数据

def read_pairs(path):
    ts, rows = [], []
    with open(path, encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split(",")
            if len(parts) < 4:
                continue
            ts.append(int(parts[0]))
            rows.append((float(parts[1]), float(parts[2]), float(parts[3])))
    return ts, rows


def read_marks(path):
    txt = open(path, encoding="utf-8").read()
    out = {}
    for tag in ("DOWN", "UP"):
        m = re.search(tag + r"\s+event=(\d+) ms\s+recv_uptime=(\d+) ms\s+recv_elapsed=(\d+) ns", txt)
        if m:
            out[tag] = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
    # App 自报的结果
    m = re.search(r"★ 距离 ([\d.]+) cm（裸积分 ([\d.]+) cm）", txt)
    out["app_cm"] = float(m.group(1)) if m else None
    m = re.search(r"窗口 #(\d+)→#(\d+)（按住 ([\d.]+) s / 积分 ([\d.]+) s / 顺延 (\d+) ms）", txt)
    out["app_win"] = (int(m.group(1)), int(m.group(2))) if m else None
    out["app_raw"] = float(m.group(2)) if False else None
    return out


def locate(ts, target):
    lo, hi = 0, len(ts) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if ts[mid] < target:
            lo = mid + 1
        else:
            hi = mid
    return lo


# ----------------------------------------------------------------- 体检

def rolling_std(rows, base):
    """逐点局部标准差（0.2 s 窗，逐轴取最大）。返回同长度列表。"""
    n = len(rows)
    w = max(3, int(0.2 / base))
    out = [0.0] * n
    half = w // 2
    for i in range(n):
        a = max(0, i - half)
        b = min(n, i + half + 1)
        m = [0.0, 0.0, 0.0]
        for k in range(a, b):
            for x in range(3):
                m[x] += rows[k][x]
        c = b - a
        s = [0.0, 0.0, 0.0]
        for k in range(a, b):
            for x in range(3):
                d = rows[k][x] - m[x] / c
                s[x] += d * d
        out[i] = max(math.sqrt(s[x] / c) for x in range(3))
    return out


def frac_static(rows, a, b, thr):
    if b - a < 3:
        return 1.0
    sd = rolling_std(rows[a:b], 0.009514)
    return sum(1 for v in sd if v <= thr) / len(sd)


def report(d, verbose=True):
    acc_ts, acc = read_pairs(os.path.join(d, "acc.csv"))
    gyr_ts, gyr = read_pairs(os.path.join(d, "gyr.csv"))
    marks = read_marks(os.path.join(d, "meta.txt"))
    name = os.path.basename(d)
    n = min(len(acc), len(gyr))
    acc, gyr = acc[:n], gyr[:n]
    dts = [acc_ts[i + 1] - acc_ts[i] for i in range(min(len(acc_ts), n) - 1)]
    dt = statistics.median(dts) / 1e9
    dn = marks.get("DOWN", (None, None, None))[2]
    up = marks.get("UP", (None, None, None))[2]
    if dn is None or up is None:
        print(f"{name}: 无触摸标记")
        return None
    lo, hi = locate(acc_ts, dn), locate(acc_ts, up)

    # 全记录里最安静的一段 = 本机真实噪声地板
    sd_all = rolling_std(acc, dt)
    floor = sorted(v for v in sd_all if v > 0)[max(0, len(sd_all) // 20)]
    static_thr = max(4.0 * floor, 0.02)

    # 现役校准段：DOWN 之前 2 s
    cal_n = max(20, min(int(2.0 / dt), lo))
    cal_lo = max(0, lo - cal_n)
    gyro_b = tuple(statistics.fmean(gyr[i][a] for i in range(cal_lo, lo)) for a in range(3))
    a_mean = tuple(statistics.fmean(acc[i][a] for i in range(cal_lo, lo)) for a in range(3))
    cal_sigma = max(statistics.pstdev([acc[i][a] for i in range(cal_lo, lo)]) for a in range(3))
    cal_stat = frac_static(acc, cal_lo, lo, static_thr)

    prof = P.DeviceProfile(dt, gyro_b, vunit(a_mean), vnorm(a_mean), cal_sigma,
                           max(statistics.pstdev([gyr[i][a] for i in range(cal_lo, lo)]) for a in range(3)),
                           cal_n)
    p = P.Params(static_mode="adaptive", assumed_dt=dt)
    truth = {"N": n - 1, "acc": acc, "gyr": gyr, "dt": dt, "n_pre": lo, "n_mot": hi - lo}
    qs, _ = P.solve_attitude(truth, prof.gyro_bias, prof.gravity_ref, p, dt=dt)
    a_lin = P.remove_gravity(truth, qs, g_mag=prof.gravity_mag)

    bias = math.sqrt(sum(statistics.fmean(a_lin[i][a] for i in range(cal_lo, lo)) ** 2 for a in range(3)))

    T = (hi - lo) * dt
    win_stat = frac_static(acc, lo, hi, static_thr)
    # 窗口内 |a_lin| 与陀螺
    al_mag = [vnorm(a_lin[i]) for i in range(lo, hi + 1)]
    gy_mag = [vnorm(gyr[i]) for i in range(lo, hi + 1)]
    al_q90 = sorted(al_mag)[int(len(al_mag) * 0.9)]
    gy_q90 = sorted(gy_mag)[int(len(gy_mag) * 0.9)]

    v, pos = P.integrate_vec(a_lin, lo, hi, dt=dt)
    bare = vnorm(pos[hi - lo])
    # 强制去趋势
    dr = [v[hi - lo][x] / T for x in range(3)]
    a2 = [tuple(a_lin[i][x] - dr[x] for x in range(3)) for i in range(lo, hi + 1)]
    _, p2 = P.integrate_vec(a2, 0, hi - lo, dt=dt)
    det = vnorm(p2[hi - lo])

    d0, ve, used, thr = P.integrate_cond_zupt3d(a_lin, lo, hi, dt, bias_bound=bias,
                                                sigma=cal_sigma, k_sigma=5.0)
    print(f"=== {name} ===")
    print(f"  按住 {T:.2f} s  窗口 #{lo}→#{hi}  n={n}  App 自报 "
          f"{marks['app_cm']} cm  （窗口 {marks['app_win']}）")
    print(f"  校准段 #{cal_lo}→#{lo}  |g|={prof.gravity_mag:.4f}（真值≈9.91）  σ_a={cal_sigma:.4f}"
          f"  |b|={bias:.4f}  陀螺零偏={tuple(round(x,6) for x in gyro_b)}")
    print(f"  ★ 校准段里真正静止的比例 = {cal_stat*100:5.1f}%   本机噪声地板={floor:.5f}  判阈={static_thr:.5f}")
    print(f"  ★ 窗口内真正静止的比例   = {win_stat*100:5.1f}%   |a_lin| 90分位={al_q90:.4f}  |gyr| 90分位={gy_q90:.5f}")
    print(f"  v_end={vnorm(v[hi-lo]):.4f}  thr={thr:.4f}  裸积分={bare*100:8.2f} cm  "
          f"强制去趋势={det*100:8.2f} cm  现役(条件)={d0*100:8.2f} cm  {'去趋势' if used else '裸积分'}")
    if verbose:
        # 静/动态时间线（每 0.25 s 一格）
        step = max(1, int(0.25 / dt))
        line = []
        for i in range(0, n, step):
            j = min(n, i + step)
            sd = rolling_std(acc[i:j], dt)
            line.append("." if (sd and max(sd) <= static_thr) else "#")
        print("  时间线（. = 静止, # = 在动）  0s → %.1fs" % ((n - 1) * dt))
        print("    " + "".join(line)[:150])
        dn_mark = max(0, lo // step)
        up_mark = max(0, hi // step)
        pfx = " " * 4
        print("    " + pfx + " " * dn_mark + "D" + " " * max(0, up_mark - dn_mark - 1) + "U")
    print()
    return dict(name=name, T=T, cal_stat=cal_stat, win_stat=win_stat,
                g=prof.gravity_mag, sigma=cal_sigma, b=bias, bare=bare * 100,
                det=det * 100, cond=d0 * 100, used=used, ve=vnorm(v[hi - lo]), thr=thr,
                app=marks["app_cm"])


def main():
    args = sys.argv[1:]
    if not args or args[0] == "--all":
        dirs = sorted(os.path.join(CAP_ROOT, x) for x in os.listdir(CAP_ROOT))
    else:
        dirs = args
    rows = []
    for d in dirs:
        r = report(d)
        if r:
            rows.append(r)
    print("=" * 118)
    print(f"{'采集':>16}{'T(s)':>7}{'校准静%':>9}{'窗口静%':>9}{'|g|':>8}{'σ_a':>8}{'|b|':>8}"
          f"{'裸积分':>10}{'去趋势':>10}{'现役':>10}{'App':>9}")
    for r in rows:
        print(f"{r['name']:>16}{r['T']:7.2f}{r['cal_stat']*100:9.1f}{r['win_stat']*100:9.1f}"
              f"{r['g']:8.4f}{r['sigma']:8.4f}{r['b']:8.4f}{r['bare']:10.2f}{r['det']:10.2f}"
              f"{r['cond']:10.2f}{(r['app'] if r['app'] is not None else float('nan')):9.2f}")


if __name__ == "__main__":
    main()
