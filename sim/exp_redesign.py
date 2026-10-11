"""重新推理：姿态参考系 + 窗口裁剪 + 无条件 ZUPT 三项改动的效果。

用户反馈：「按着不动都能测出来很长一段距离」「对姿态要求也很严格」。

用法: python exp_redesign.py
"""
import math
import os
import statistics

import pipeline as P
from pipeline import (G, q_delta, q_from_two_vectors, qnorm, qrotT, qmul,
                      vcross, vnorm, vscale, vsub, vunit)
import diag as D

ALPHA = 0.05
GATE_RAD = 0.3
MARGIN_S = 0.2


# ---------------------------------------------------------------- 工具

def rolling_std_max(rows, dt, win_s=0.2):
    n = len(rows)
    w = max(3, int(win_s / dt))
    out = [0.0] * n
    h = w // 2
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        c = b - a
        best = 0.0
        for x in range(3):
            m = sum(rows[k][x] for k in range(a, b)) / c
            best = max(best, math.sqrt(sum((rows[k][x] - m) ** 2 for k in range(a, b)) / c))
        out[i] = best
    return out


def static_mask(acc, gyr, gb, dt):
    sd = rolling_std_max(acc, dt)
    floor = sorted(v for v in sd if v > 0)[max(0, len(sd) // 20)]
    thr = max(4.0 * floor, 0.02)
    om = [vnorm(vsub(gyr[i], gb)) for i in range(len(gyr))]
    ofloor = sorted(om)[max(0, len(om) // 20)]
    othr = max(4.0 * ofloor, 0.01)
    w = max(3, int(0.2 / dt))
    flags = []
    for i in range(len(acc)):
        a, b = max(0, i - w + 1), i + 1
        flags.append(max(sd[a:b]) <= thr and max(om[a:b]) <= othr)
    return flags, floor, thr, ofloor, othr


def solve(acc, gyr, gb, g_ref_body, flags, dt, ref):
    """ref='body' = 现役（锁死在校准姿态的体坐标参考）；
       ref='measured' = 修正（用当前实测加速度方向，且只在静止处修正）。"""
    N = len(acc) - 1
    q0 = q_from_two_vectors(g_ref_body, (0.0, 0.0, 1.0))
    q = q0
    qs = [q]
    skipped = 0
    for i in range(N):
        w = vsub(gyr[i], gb)
        q_pred = qnorm(qmul(q, q_delta(w, dt)))
        g_pred_body = qrotT(q_pred, (0.0, 0.0, 1.0))
        if ref == "measured":
            if not flags[i]:
                q = q_pred
                skipped += 1
                qs.append(q)
                continue
            g_ref_i = vunit(acc[i])
        else:
            g_ref_i = g_ref_body
        e = vcross(g_pred_body, g_ref_i)
        if vnorm(e) > GATE_RAD:
            q = q_pred
            skipped += 1
        else:
            q = qnorm(qmul(q_pred, q_delta(vscale(e, -ALPHA), dt)))
        qs.append(q)
    return qs, skipped


def a_lin_of(acc, qs, g_mag):
    return [P.qrot(qs[i], acc[i])[:2] + (P.qrot(qs[i], acc[i])[2] - g_mag,) for i in range(len(acc))]


def detrend_dist(a_lin, lo, hi, dt):
    T = (hi - lo) * dt
    v = [0.0, 0.0, 0.0]
    for i in range(lo, hi):
        for x in range(3):
            v[x] += 0.5 * (a_lin[i][x] + a_lin[i + 1][x]) * dt
    dr = [v[x] / T if T > 1e-9 else 0.0 for x in range(3)]
    v = [0.0, 0.0, 0.0]
    p = [0.0, 0.0, 0.0]
    for i in range(lo, hi):
        for x in range(3):
            aa = 0.5 * ((a_lin[i][x] - dr[x]) + (a_lin[i + 1][x] - dr[x]))
            v2 = v[x] + aa * dt
            p[x] += 0.5 * (v[x] + v2) * dt
            v[x] = v2
    return math.sqrt(sum(x * x for x in p)), vnorm(v)


def bare_dist(a_lin, lo, hi, dt):
    v = [0.0, 0.0, 0.0]
    p = [0.0, 0.0, 0.0]
    for i in range(lo, hi):
        for x in range(3):
            aa = 0.5 * (a_lin[i][x] + a_lin[i + 1][x])
            v2 = v[x] + aa * dt
            p[x] += 0.5 * (v[x] + v2) * dt
            v[x] = v2
    return math.sqrt(sum(x * x for x in p))


# ---------------------------------------------------------------- 主流程

def run(name, truth_cm=None, verbose=True):
    d = os.path.join(D.CAP_ROOT, name)
    acc_ts, acc = D.read_pairs(os.path.join(d, "acc.csv"))
    gyr_ts, gyr = D.read_pairs(os.path.join(d, "gyr.csv"))
    marks = D.read_marks(os.path.join(d, "meta.txt"))
    if marks.get("DOWN") is None or marks.get("UP") is None:
        return None
    n = min(len(acc), len(gyr))
    acc, gyr = acc[:n], gyr[:n]
    dt = statistics.median([acc_ts[i + 1] - acc_ts[i] for i in range(n - 1)]) / 1e9
    lo, hi = D.locate(acc_ts, marks["DOWN"][2]), D.locate(acc_ts, marks["UP"][2])
    flags, floor, thr, ofloor, othr = static_mask(acc, gyr, tuple(
        statistics.fmean(gyr[i][a] for i in range(max(0, lo - 200), lo)) for a in range(3)), dt)

    # 校准段：在 DOWN 之前 4 s 内取最长的一整段静止
    a4 = max(0, lo - int(4.0 / dt))
    runs = []
    i = a4
    while i < lo:
        if flags[i]:
            j = i
            while j + 1 < lo and flags[j + 1]:
                j += 1
            runs.append((i, j))
            i = j + 1
        else:
            i += 1
    if runs:
        runs.sort(key=lambda r: r[1] - r[0])
        cal_lo, cal_hi = runs[-1][0], runs[-1][1] + 1
    else:
        cal_lo, cal_hi = max(0, lo - int(2.0 / dt)), lo
    cn = cal_hi - cal_lo
    gb = tuple(statistics.fmean(gyr[i][a] for i in range(cal_lo, cal_hi)) for a in range(3))
    am = tuple(statistics.fmean(acc[i][a] for i in range(cal_lo, cal_hi)) for a in range(3))
    g_mag = vnorm(am)
    g_ref = vunit(am)
    sig = max(statistics.pstdev([acc[i][a] for i in range(cal_lo, cal_hi)]) for a in range(3))

    # 运动段（按钮窗内第一/最后一个"在动"的样本）
    mot = [i for i in range(lo, hi + 1) if not flags[i]]
    if mot:
        m0, m1 = mot[0], mot[-1]
        margin = int(MARGIN_S / dt)
        w0, w1 = max(lo, m0 - margin), min(hi, m1 + margin)
        # 端点必须落在静止处
        while w0 < m0 and not flags[w0]:
            w0 += 1
        while w1 > m1 and not flags[w1]:
            w1 -= 1
    else:
        w0, w1 = None, None

    out = {"name": name, "truth": truth_cm, "lo": lo, "hi": hi, "dt": dt}
    if verbose:
        print(f"=== {name} ===  真值 {truth_cm} cm")
        print(f"  噪声地板 acc {floor:.5f} 陀螺 {ofloor:.5f}  校准段 #{cal_lo}→#{cal_hi}"
              f"（{cn*dt:.2f} s，静止 {sum(1 for i in range(cal_lo,cal_hi) if flags[i])/cn*100:.0f}%）"
              f"  |g|={g_mag:.4f}  σ_a={sig:.4f}")
        if mot:
            print(f"  按钮窗 #{lo}→#{hi}（{(hi-lo)*dt:.2f} s）内运动段 #{m0}→#{m1}"
                  f"（{(m1-m0)*dt:.2f} s）→ 裁剪窗 #{w0}→#{w1}（{(w1-w0)*dt:.2f} s）")
        else:
            print(f"  按钮窗 #{lo}→#{hi}（{(hi-lo)*dt:.2f} s）内**没有任何运动** → 距离应为 0")

    for ref in ("body", "measured"):
        qs, skip = solve(acc, gyr, gb, g_ref, flags, dt, ref)
        al = a_lin_of(acc, qs, g_mag)
        # 校准段残余零偏（矢量模）
        bb = vnorm(tuple(statistics.fmean(al[i][x] for i in range(cal_lo, cal_hi)) for x in range(3)))
        btn_det, _ = detrend_dist(al, lo, hi, dt)
        btn_bare = bare_dist(al, lo, hi, dt)
        if w0 is not None:
            trim_det, _ = detrend_dist(al, w0, w1, dt)
        else:
            trim_det = 0.0
        out[ref] = dict(skip=skip, bias=bb, btn_bare=btn_bare * 100, btn_det=btn_det * 100,
                        trim=trim_det * 100)
        if verbose:
            t = f"  {ref:>8}"
            print(f"{t}  跳过加速度修正 {skip:4d}  |b|={bb:7.4f}   "
                  f"按钮窗:裸 {btn_bare*100:8.2f} / 去趋势 {btn_det*100:8.2f}   "
                  f"裁剪窗去趋势 {trim_det*100:8.2f} cm")
    if verbose:
        print()
    return out


def main():
    names = sorted(x for x in os.listdir(D.CAP_ROOT) if os.path.isdir(os.path.join(D.CAP_ROOT, x)))
    rows = []
    for nm in names:
        r = run(nm, verbose=False)
        if r is not None:
            rows.append(r)
    print("=" * 132)
    print("现役(body 参考) vs 修正(measured 参考)；裸积分 vs 无条件去趋势 vs 裁剪窗去趋势（cm）")
    print(f"{'采集':>16}{'按钮窗s':>9}{'|b|现役':>9}{'现役裸':>9}{'现役去趋势':>11}"
          f"{'|b|修正':>9}{'修正裸':>9}{'修正去趋势':>11}{'修正+裁剪':>10}{'App':>9}")
    for r in rows:
        b, m = r["body"], r["measured"]
        print(f"{r['name']:>16}{(r['hi']-r['lo'])*r['dt']:9.2f}{b['bias']:9.4f}{b['btn_bare']:9.2f}"
              f"{b['btn_det']:11.2f}{m['bias']:9.4f}{m['btn_bare']:9.2f}{m['btn_det']:11.2f}"
              f"{m['trim']:10.2f}")


if __name__ == "__main__":
    main()
