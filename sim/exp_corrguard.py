"""「去趋势修正量过大 ⟹ 该条不可信」—— 一个候选守卫的检验。

背景（2026-10-17）：26 条手操 20 cm 采集里，误差与**运动时长**分开得很干净：

    ≤3.77 s 的 12 条   误差全在 ±0.9 cm 内（= 尺子自身精度）
    第一个跑出 ±1 cm 的是 3.85 s
    >3.8 s 的 14 条    11 条跑出，最大 +7.81 / −6.13

而 sqrt 里那个"去趋势修正量"（bare − engine）在长窗上会变得很大。于是有个自然的
假设：**修正量本身大，就是该条不可信的信号**——因为两个锚点只给 2 个约束，
线性残差正好 2 个参数，**没有冗余、没有自检能力**；修正量越大，说明模型在
越用力地跟数据搏斗。

本脚本把这个假设变成数字，并检查一个简单守卫的收益/代价：

    「若 |修正量| 超过阈值 → 不给数、判为不可信」

⚠ 两条必须写在前面的话：

1. **这是同批内（in-sample）评估。** 阈值是在这 26 条上挑的，必然偏乐观。
   把它接进引擎之前，**必须先在一批完全没有参与定阈值的采集上复核**。
   本项目的翻车史（见 sim/README.md §0.5）里有三次都是"在一小批上看着成立
   就升级成结论"。这里刻意**只做诊断候选，不动引擎**。
2. 阈值优先用**比例**（|修正|/bare）而不是绝对值 —— 绝对阈值会随距离漂移，
   20 cm 上合适的数到了 40 cm 上就不合适。比例阈值更可能跨距离成立。

用法：
    python sim/exp_corrguard.py
"""

import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2
from exp_real import load, ROOT, HAND20, TRUTH_CM


def collect():
    rows = []
    for nm in sorted(HAND20):
        acc, gyr, dn, up, dt = load(os.path.join(ROOT, nm))
        m = measure2.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                             head_r=None, gate_motion=True)
        if not m["ok"]:
            continue
        truth = TRUTH_CM[nm]
        bare = m["bareM"] * 100
        eng = m["distanceM"] * 100
        corr = abs(bare - eng)              # 去趋势实际扣掉的量
        rows.append({
            "nm": nm, "movS": m["motionS"], "seg": m["nMerged"],
            "bare": bare, "eng": eng, "err": eng - truth,
            "absErr": abs(eng - truth),
            "corr": corr, "corrPct": corr / bare * 100.0,
        })
    return rows


def pearson(a, b):
    ma, mb = statistics.fmean(a), statistics.fmean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5
    return num / den if den else 0.0


def main():
    rows = collect()
    n = len(rows)
    print("真值 20 cm（手操，自报 ±1 cm）  有效 %d 条" % n)
    print()
    print("%-18s %5s %3s %7s %7s %7s %8s %8s"
          % ("采集", "运动s", "段", "裸积", "引擎", "误差", "修正cm", "修正%"))
    print("-" * 76)
    for r in sorted(rows, key=lambda r: r["corr"]):
        mark = "  <<<" if r["absErr"] > 1.0 else ""
        print("%-18s %5.2f %3d %7.2f %7.2f %+7.2f %8.2f %7.1f%%%s"
              % (r["nm"], r["movS"], r["seg"], r["bare"], r["eng"],
                 r["err"], r["corr"], r["corrPct"], mark))
    print("-" * 76)

    good = [r for r in rows if r["absErr"] <= 1.0]
    bad = [r for r in rows if r["absErr"] > 1.0]
    print()
    print("好区（|误差| ≤ 1 cm） n=%d   修正量中位 %.2f cm (%.1f%%)  最大 %.2f cm (%.1f%%)"
          % (len(good), statistics.median([r["corr"] for r in good]),
             statistics.median([r["corrPct"] for r in good]),
             max(r["corr"] for r in good), max(r["corrPct"] for r in good)))
    print("尾部（|误差| > 1 cm） n=%d   修正量中位 %.2f cm (%.1f%%)  最大 %.2f cm (%.1f%%)"
          % (len(bad), statistics.median([r["corr"] for r in bad]),
             statistics.median([r["corrPct"] for r in bad]),
             max(r["corr"] for r in bad), max(r["corrPct"] for r in bad)))
    print()
    print("corr(|误差|, 修正量)     = %+.3f" % pearson([r["absErr"] for r in rows], [r["corr"] for r in rows]))
    print("corr(|误差|, 运动时长)   = %+.3f" % pearson([r["absErr"] for r in rows], [r["movS"] for r in rows]))

    print()
    print("守卫「修正量超阈值 → 判不可信」的收益 / 代价（in-sample，仅供定候选）：")
    print("  比例阈值：拦下尾部 / 误伤好区")
    for th in (10, 15, 20, 25, 30):
        caught = sum(1 for r in bad if r["corrPct"] > th)
        hurt = sum(1 for r in good if r["corrPct"] > th)
        print("    修正 > %2d%% :  尾部 %2d/%d 拦下，好区 %d/%d 误伤"
              % (th, caught, len(bad), hurt, len(good)))
    print("  绝对阈值：")
    for th in (1.0, 2.0, 3.0, 4.0):
        caught = sum(1 for r in bad if r["corr"] > th)
        hurt = sum(1 for r in good if r["corr"] > th)
        print("    修正 > %.1fcm :  尾部 %2d/%d 拦下，好区 %d/%d 误伤"
              % (th, caught, len(bad), hurt, len(good)))

    print()
    print("⚠ 同批内挑的阈值必然偏乐观。接引擎前必须在一批全新采集上复核；")
    print("  在那之前，它只配当**诊断行**显示，不参与给数/拒答。")


if __name__ == "__main__":
    main()
