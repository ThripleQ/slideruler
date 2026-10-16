"""实验 16：姿态滤波的两个结构性缺陷（真机硬真值 15 cm ×3）

现象（20261011_143719/143725/143734，姿态"比较随意"）：
  静止段去重力后残留 2.4~3.3 m/s² 的**常值**，相当于 ~19° 的重力泄漏。
  这么大的常值不可能是加速度计零偏（真实零偏 ~0.05 m/s²）。

根因（两条，互相掩盖）：
  F1 初始化错位：初始姿态取自校准段（样本 172~301）的重力方向，姿态积分
                 却从样本 0 开始。若手机在校准段之前转过（实测 21°），
                 初始姿态对样本 0 就是错的，且这个错一直带进窗口。
  F2 门控死锁：  |e| > GATE_RAD(0.30 rad ≈ 17°) 时跳过修正。误差一旦超过
                 门限就再也拉不回来 —— 它保护的正是它自己没能防住的东西。

本脚本不改动 measure2.py，用三个变体对照：
  base    现状
  fix1    只修 F1（从校准段起点开始积分姿态）
  fix12   同时修 F1 + F2（静止段不再用方向门限，改用 |acc| 模长检验）
"""
import argparse
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pipeline import (q_delta, q_from_two_vectors, qmul, qnorm, qrot, qrotT,
                      vcross, vnorm, vscale, vsub, vunit)
import measure2 as M
from exp_real import load, ROOT

# 硬真值：用户用 15 cm 直尺，2026-10-11 14:37 的最后三次
TRUTH = {"20261011_143719": 15.0, "20261011_143725": 15.0,
         "20261011_143734": 15.0}


def attitude(acc, gyr, gb, g_ref_body, flags, dt, i0=0, alpha=M.ALPHA,
             gate=M.GATE_RAD, gate_mode="dir", g_mag=None):
    """solve_attitude_v2 的可配置版。

    i0        从第 i0 个样本开始积分（i0 之前填初始姿态）
    gate_mode "dir" 方向门限（现状，会死锁）
              "mag" |acc| 模长检验：模长偏离 |g| 太多才不修正（静止段 acc 就是重力）
              "off" 静止段一律修正
    """
    q = q_from_two_vectors(g_ref_body, (0.0, 0.0, 1.0))
    qs = [q] * i0
    tol = 0.25 * (g_mag or 9.8)
    for i in range(i0, len(acc) - 1):
        q_pred = qnorm(qmul(q, q_delta(vsub(gyr[i], gb), dt)))
        if flags[i]:
            e = vcross(qrotT(q_pred, (0.0, 0.0, 1.0)), vunit(acc[i]))
            if gate_mode == "dir":
                ok = vnorm(e) <= gate
            elif gate_mode == "mag":
                ok = abs(vnorm(acc[i]) - (g_mag or 9.8)) <= tol
            else:
                ok = True
            q = (q_pred if not ok
                 else qnorm(qmul(q_pred, q_delta(vscale(e, -alpha), dt))))
        else:
            q = q_pred
        qs.append(q)
    while len(qs) < len(acc):
        qs.append(q)
    return qs


def run(nm, fix1=True, fix2=False, gate_mode="mag"):
    acc, gyr, dn, up, dt = load(os.path.join(ROOT, nm))
    m = M.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                  head_r=None, gate_motion=False)
    # ⚠ 不能用 `m["ok"]`：模型敏感现在直接拒答（measure2.py §8），而 143734
    # 正是那条 —— 本实验要诊断它，跳过就没意义了。拒答返回里 cal/run0/run1/flags/gMag
    # 都在，够这里继续算；只有"找不到两端锚点"的早期拒答才真算不了。
    if "run0" not in m:
        return None, m["reason"], None
    ca, cb = m["cal"]
    r0, r1 = m["run0"], m["run1"]
    gb2 = tuple(statistics.fmean(gyr[i][a] for i in range(ca, cb)) for a in range(3))
    am = tuple(statistics.fmean(acc[i][a] for i in range(ca, cb)) for a in range(3))
    g_mag = m["gMag"]
    qs = attitude(acc, gyr, gb2, vunit(am), m["flags"], dt,
                  i0=(ca if fix1 else 0),
                  gate_mode=("dir" if not fix2 else gate_mode), g_mag=g_mag)
    a_lin = M.to_linear(acc, qs, g_mag)
    dist, vc, p, drift, spread = M.zupt_two_anchor(
        a_lin, m["window"][0], m["window"][1], r0, r1, dt)
    # 静止段的残余（修复后应当 ≈ 0）
    res0 = vnorm([statistics.fmean(a_lin[i][k] for i in range(*r0)) for k in range(3)])
    res1 = vnorm([statistics.fmean(a_lin[i][k] for i in range(*r1)) for k in range(3)])
    return dist * 100, None, (res0, res1, drift)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", default=",".join(TRUTH))
    a = ap.parse_args()
    names = a.names.split(",")

    print("%-16s %6s | %9s %9s %9s" % ("采集", "真值", "现状", "只修F1", "F1+F2"))
    print("-" * 60)
    rows = []
    for nm in names:
        t = TRUTH.get(nm)
        out = []
        for fx1, fx2 in [(False, False), (True, False), (True, True)]:
            d, why, extra = run(nm, fix1=fx1, fix2=fx2)
            out.append((d, why, extra))
        base = out[0][0]
        f1 = out[1][0]
        f12 = out[2][0]

        # 拒答（143734 现在由模型敏感拒答）要显示"拒答"而不是 0.00 ——
        # 否则那一列看起来像"算法把它测没了"（−100%），实际是引擎明确拒绝给数。
        def cell(v):
            return "     拒答" if v is None else "%8.2f" % v

        print("%-16s %6.1f | %s %s %s" % (nm, t or -1, cell(base), cell(f1), cell(f12)))
        for tag, (d, why, extra) in zip(("现状", "F1", "F1+F2"), out):
            if extra:
                print("      %-6s 静止段残余 %.3f / %.3f m/s²  锚点漂移 %.3f m/s"
                      % (tag, extra[0], extra[1], extra[2]))
        rows.append((nm, t, base, f1, f12))

    print("-" * 60)
    # 统计只涵盖三条都给出数的那些（拒答的不能按 0 算）。
    got = [r for r in rows if r[1] and all(r[i] is not None for i in (2, 3, 4))]
    if got:
        for i, tag in enumerate(("现状", "只修F1", "F1+F2"), start=2):
            errs = [(r[i] - r[1]) / r[1] * 100 for r in got]
            print("%-6s 相对误差 %s   平均 %+.2f%%  最大 |%.2f%%|"
                  % (tag, "  ".join("%+.1f%%" % e for e in errs),
                     statistics.fmean(errs), max(abs(e) for e in errs)))


if __name__ == "__main__":
    main()
