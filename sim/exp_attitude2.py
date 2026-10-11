"""实验 18：姿态双锚点 —— 姿态和速度必须由不同的东西定标

死结在于把两件事混为一谈了：

  姿态（重力方向）  静止段 acc = 重力，**匀速段 acc 也 = 重力**（比力 = 重力，
                    a_true = 0，伽利略）。所以匀速段一样能锚定姿态。
  速度/位置         只有**真静止**才有 v = 0。匀速段 v ≠ 0，把它强行归零
                    就会把真实位移扣掉（143734 因此被扣掉 14.8 cm，答案剩 1.59）。

现状（双锚点 ZUPT）用同一段同时定这两样，于是：
  · 终点在匀速 → 速度锚点无效 → 位移被吞
  · 运动段姿态漂移（陀螺零偏游走）无人管 → 泄漏 ∝ ½·g·b_g·T³

本实验把它们拆开：
  1. 起点静止段定初始姿态与加速度零偏 b_a
  2. 终点段（静止或匀速）的 acc 方向反解**陀螺零偏残量** b_g
     —— 把运动段的姿态漂移当成常值角速度误差，两端一拉就平直
  3. 用修正后的姿态重算 a_lin，只扣 b_a，**不强行把末端速度归零**

判定：真值 15 cm ×3，谁最接近。
"""
import argparse
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pipeline import (q_delta, q_from_two_vectors, qmul, qnorm, qrot, qrotT,
                      qconj, vcross, vnorm, vscale, vsub, vunit)
import measure2 as M
from exp_real import load, ROOT
from exp_attitude import TRUTH


def integrate_q(q0, gyr, gb, i0, i1, dt):
    """陀螺积分，返回 qs[i0..i1]。"""
    q = q0
    qs = [q]
    for i in range(i0, i1):
        q = qnorm(qmul(q, q_delta(vsub(gyr[i], gb), dt)))
        qs.append(q)
    return qs


def rotvec(q):
    """四元数 → 旋转向量（弧度，方向 = 旋转轴）。"""
    q = qnorm(q)
    w = max(-1.0, min(1.0, q[0]))
    ang = 2.0 * math.atan2(vnorm(q[1:]), abs(q[0]))
    if ang < 1e-12:
        return (0.0, 0.0, 0.0)
    if q[0] < 0:
        s = -ang / vnorm(q[1:])
    else:
        s = ang / vnorm(q[1:])
    return tuple(q[k] * s for k in (1, 2, 3))


def solve_bg(q_start, gyr, gb, i0, i1, dt, g_end_body):
    """反解陀螺零偏残量 b_g，使积分到 i1 时预测重力方向 = g_end_body。

    预测重力方向（体系）由姿态决定；实测的 g_end_body 来自终点段 acc 均值。
    两者之差是一个小旋转，把它均摊到 (i1-i0)·dt 上就是常值角速度误差。
    """
    qs = integrate_q(q_start, gyr, gb, i0, i1, dt)
    q_pred = qs[-1]
    a_pred = qrotT(q_pred, (0.0, 0.0, 1.0))       # 预测重力方向（体系）
    e = vcross(a_pred, g_end_body)                 # 把预测转向实测的旋转轴
    s = vnorm(e)
    if s < 1e-12:
        return (0.0, 0.0, 0.0), 0.0
    ang = math.atan2(s, sum(a_pred[k] * g_end_body[k] for k in range(3)))
    # 体系内的旋转向量
    rv = tuple(e[k] / s * ang for k in range(3))
    T = (i1 - i0) * dt
    return tuple(rv[k] / T for k in range(3)), math.degrees(ang)


def run(nm, mode="win"):
    acc, gyr, dn, up, dt = load(os.path.join(ROOT, nm))
    m = M.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                  head_r=None, gate_motion=False)
    if not m["ok"]:
        return None
    ca, cb = m["cal"]
    r0, r1 = m["run0"], m["run1"]
    w0, w1 = m["window"]
    mo0, mo1 = m["motion"]

    gb2 = tuple(statistics.fmean(gyr[i][a] for i in range(ca, cb)) for a in range(3))
    am0 = tuple(statistics.fmean(acc[i][a] for i in range(ca, cb)) for a in range(3))
    am1 = tuple(statistics.fmean(acc[i][a] for i in range(*r1)) for a in range(3))
    g_mag = m["gMag"]
    q0 = q_from_two_vectors(vunit(am0), (0.0, 0.0, 1.0))

    # 反解陀螺零偏残量（在体系内）
    t_end = (r1[0] + r1[1]) // 2
    bg, ang_deg = solve_bg(q0, gyr, gb2, ca, t_end, dt, vunit(am1))
    gb_all = tuple(gb2[k] + bg[k] for k in range(3))

    # 用修正后的陀螺重算姿态（静止段仍然用 acc 修正，运动段纯外推）
    q = q0
    qs = [q0] * ca
    for i in range(ca, len(acc) - 1):
        q_pred = qnorm(qmul(q, q_delta(vsub(gyr[i], gb_all), dt)))
        if m["flags"][i] and abs(vnorm(acc[i]) - g_mag) <= 0.25 * g_mag:
            e = vcross(qrotT(q_pred, (0.0, 0.0, 1.0)), vunit(acc[i]))
            q = qnorm(qmul(q_pred, q_delta(vscale(e, -M.ALPHA), dt)))
        else:
            q = q_pred
        qs.append(q)
    while len(qs) < len(acc):
        qs.append(qs[-1])

    a = M.to_linear(acc, qs, g_mag)
    b_a = tuple(statistics.fmean(a[i][k] for i in range(*r0)) for k in range(3))

    def integ(i0, i1):
        v = [0.0, 0.0, 0.0]
        p = [0.0, 0.0, 0.0]
        for i in range(i0 + 1, i1 + 1):
            vv = [v[k] + 0.5 * ((a[i - 1][k] - b_a[k]) + (a[i][k] - b_a[k])) * dt
                  for k in range(3)]
            p = [p[k] + 0.5 * (v[k] + vv[k]) * dt for k in range(3)]
            v = vv
        return vnorm(v), vnorm(p)

    v_end, p_win = integ(w0, w1)
    v_end_m, p_mot = integ(mo0, mo1)

    # 对照：现状双锚点（旧姿态）
    return dict(旧双锚点=m["distanceM"] * 100,
                姿态锚_窗口=p_win * 100,
                姿态锚_运动段=p_mot * 100,
                _diag=dict(陀螺残量rad_s=vnorm(bg), 终点姿态差=ang_deg,
                           末端速度=v_end_m, 运动T=(mo1 - mo0) * dt,
                           零偏b_a=vnorm(b_a),
                           终点残余=vnorm(tuple(
                               statistics.fmean(a[i][k] for i in range(*r1))
                               for k in range(3)))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", default=",".join(TRUTH))
    a = ap.parse_args()
    names = a.names.split(",")
    keys = ["旧双锚点", "姿态锚_窗口", "姿态锚_运动段"]
    rows = []
    print("%-16s %5s | %s" % ("采集", "真值", "  ".join("%13s" % k for k in keys)))
    print("-" * (24 + 15 * len(keys)))
    for nm in names:
        r = run(nm)
        if r is None:
            print(nm, "拒答"); continue
        rows.append((nm, TRUTH.get(nm), r))
        print("%-16s %5.1f | %s" % (nm, TRUTH.get(nm) or -1,
              "  ".join("%13.2f" % r[k] for k in keys)))
        d = r["_diag"]
        print("      %-14s 陀螺残量 %.4f rad/s  终点姿态差 %.2f°  末端v %.3f  "
              "运动%.2fs  |b_a|=%.3f  终点残余 %.3f"
              % ("", d["陀螺残量rad_s"], d["终点姿态差"], d["末端速度"],
                 d["运动T"], d["零偏b_a"], d["终点残余"]))
    print("-" * (24 + 15 * len(keys)))
    for k in keys:
        errs = [(r[2][k] - r[1]) / r[1] * 100 for r in rows if r[1]]
        if errs:
            print("%-14s 相对误差 %s  平均 %+.2f%%  最大 |%.2f%%|"
                  % (k, "  ".join("%+.1f%%" % e for e in errs),
                     statistics.fmean(errs), max(abs(e) for e in errs)))


if __name__ == "__main__":
    main()
