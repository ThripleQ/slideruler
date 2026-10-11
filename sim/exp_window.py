"""实验 17：窗口该多长？—— 硬真值 15 cm ×3 上的对照

上一发（exp_attitude.py）修掉了姿态初始化错位（静止段残余 3.3 → 0.002 m/s²），
但距离没变（11.64→11.78）。说明那个巨大的**常值**泄漏早被双锚点扣掉了，
它不是主误差。真正的矛盾在窗口：

  143734  窗口 0.97 s（真实运动只有 0.69 s）  裸积分 16.39  双锚点 1.59
  143725  窗口 1.78 s                        裸积分 21.21  双锚点 15.96
  143719  窗口 1.72 s                        裸积分  4.21  双锚点 11.71

三条的末端速度分别是 0.333 / 0.133 / 0.227 m/s，**都不为 0** —— 松手那一刻
手机还在动。双锚点的"两端速度归零"把这个真实速度当成漂移扣掉，
143734 因此被扣掉 14.8 cm（正好是答案）。

同时窗口里塞了 0.3~1.1 s 的静止段。静止段对位移零贡献，对误差却是完整加成
（误差 ∝ ½·b·T²）—— 这是 docs/boundary-design.md §3 早就写下的硬事实。

本实验对照四种窗口/定标组合：
  A 现状     窗口 = 静止→静止，双锚点（两端强行归零）
  B 短窗口   窗口 = 纯运动段，零偏只用**窗口外**的起点静止段估（不强行归零）
  C 短窗口+两端零偏线性插值（若终点段速度≈0 才用）
  D 按住窗   窗口 = 按下→松开

判定标准只有一条：真值 15 cm，谁最接近。
"""
import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pipeline import vnorm, vunit
import measure2 as M
from exp_real import load, ROOT
from exp_attitude import attitude, TRUTH


def prep(nm):
    acc, gyr, dn, up, dt = load(os.path.join(ROOT, nm))
    m = M.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                  head_r=None, gate_motion=False)
    if not m["ok"]:
        return None
    ca, cb = m["cal"]
    gb2 = tuple(statistics.fmean(gyr[i][a] for i in range(ca, cb)) for a in range(3))
    am = tuple(statistics.fmean(acc[i][a] for i in range(ca, cb)) for a in range(3))
    qs = attitude(acc, gyr, gb2, vunit(am), m["flags"], dt, i0=ca,
                  gate_mode="mag", g_mag=m["gMag"])
    a = M.to_linear(acc, qs, m["gMag"])
    return dict(acc=acc, gyr=gyr, dn=dn, up=up, dt=dt, m=m, a=a,
                cal=(ca, cb), r0=m["run0"], r1=m["run1"])


def integ(a, i0, i1, dt, b=None):
    """从 i0 积到 i1，扣除常值零偏 b。返回 (v_end, p_end)。"""
    v = [0.0, 0.0, 0.0]
    p = [0.0, 0.0, 0.0]
    bb = b or (0.0, 0.0, 0.0)
    for i in range(i0 + 1, i1 + 1):
        vv = [v[k] + 0.5 * ((a[i - 1][k] - bb[k]) + (a[i][k] - bb[k])) * dt
              for k in range(3)]
        p = [p[k] + 0.5 * (v[k] + vv[k]) * dt for k in range(3)]
        v = vv
    return vnorm(v), vnorm(p)


def bias_from(a, runs):
    idx = [i for r in runs for i in range(*r)]
    return tuple(statistics.fmean(a[i][k] for i in idx) for k in range(3))


def variants(d):
    a, dt, m = d["a"], d["dt"], d["m"]
    w0, w1 = m["window"]
    mo0, mo1 = m["motion"]
    r0, r1 = m["run0"], m["run1"]
    out = {}

    # A 现状：双锚点
    out["A 现状双锚点"] = m["distanceM"] * 100

    # B 短窗口 + 起点静止段定零偏
    b0 = bias_from(a, [r0])
    _, pB = integ(a, mo0, mo1, dt, b0)
    out["B 短窗口+起锚"] = pB * 100

    # C 短窗口 + 两端零偏线性插值（零偏随时间线性漂移）
    b1 = bias_from(a, [r1])
    t0 = statistics.fmean(range(*r0)) * dt
    t1 = statistics.fmean(range(*r1)) * dt
    tm = statistics.fmean(range(mo0, mo1)) * dt
    if abs(t1 - t0) > 1e-9:
        f = (tm - t0) / (t1 - t0)
        bC = tuple(b0[k] + f * (b1[k] - b0[k]) for k in range(3))
    else:
        bC = b0
    _, pC = integ(a, mo0, mo1, dt, bC)
    out["C 短窗口+插值"] = pC * 100

    # D 按住窗 + 起点零偏
    _, pD = integ(a, d["dn"], d["up"], dt, b0)
    out["D 按住窗"] = pD * 100

    # 诊断量
    vB, _ = integ(a, mo0, mo1, dt, b0)
    diag = dict(窗口T=m["T"], 运动T=(mo1 - mo0) * dt, 末端速度=vB,
                零偏=vnorm(b0), 末端零偏=vnorm(b1))
    return out, diag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", default=",".join(TRUTH))
    a = ap.parse_args()
    names = a.names.split(",")

    keys = None
    rows = []
    for nm in names:
        d = prep(nm)
        if d is None:
            print(nm, "拒答"); continue
        out, diag = variants(d)
        keys = list(out)
        rows.append((nm, TRUTH.get(nm), out, diag))

    print("%-16s %5s | %s" % ("采集", "真值",
          "  ".join("%14s" % k for k in keys)))
    print("-" * (24 + 16 * len(keys)))
    for nm, t, out, diag in rows:
        print("%-16s %5.1f | %s" % (nm, t or -1,
              "  ".join("%14.2f" % out[k] for k in keys)))
        print("      %-18s 窗口%.2fs 运动%.2fs 末端v%.3f |b0|=%.3f |b1|=%.3f"
              % ("", diag["窗口T"], diag["运动T"], diag["末端速度"],
                 diag["零偏"], diag["末端零偏"]))

    print("-" * (24 + 16 * len(keys)))
    for k in keys:
        errs = [(r[2][k] - r[1]) / r[1] * 100 for r in rows if r[1]]
        if errs:
            print("%-16s 相对误差 %s  平均 %+.2f%%  最大 |%.2f%%|"
                  % (k, "  ".join("%+.1f%%" % e for e in errs),
                     statistics.fmean(errs), max(abs(e) for e in errs)))


if __name__ == "__main__":
    main()
