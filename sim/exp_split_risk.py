"""逐簇求和的**证伪测试** —— 在旧语料上找它的失败模式。

`exp_split.py` 在 26 条手操 20 cm 上把 RMS 从 2.50 压到 0.44 cm。但那条数据有一个
结构上的便利：**每次"多段运动"里只有一段带位移，其余是 0.01 cm 的点按敲击**。
换一种工况，这个前提可能不成立。

## 已知的失败模式（v3 的旧伤）

静止检测分不出匀速和静止（加速度计测比力）。一次**连续滑动**在它眼里会被切成
好几段，段间的"静止"其实是匀速运动。现役 `merge_mode="press"` 取并集正是为了
挡这个（真机 122416：一次按住被切 5 段，只取一段 → 0.59 cm）。

**逐簇求和会把这些段拆开各自去趋势** —— 而每段的锚点落在"以为是静止的匀速段"
上 ⟹ 每段都扣错 ⟹ 求和可能比合并更糟。

## 本脚本做什么

在**全库**上跑「合并 vs 逐簇求和」，按批次分组看两者差多少，重点看：

1. 三条硬真值 15 cm（143719/143725/143734）—— 有真值，直接看谁更近；
2. 三条肉眼估 15 cm（100444/100459/100529）；
3. `122416` —— v3 那条"一次按住被切 5 段"的标本，预期逐簇会崩；
4. 其余采集看差值的**散布**，找出 |求和−合并| 大的那些，人工看是哪种工况。

差值大 ≠ 逐簇错（可能合并错）。所以这里只**列出来**，不下结论。

## 结果（2026-10-17）

**证伪成功 —— 逐簇求和不能上。** 全库 49 条里：

- `|求和−合并|` 中位 0.20 cm、P90 3.92、最大 **96.07**；
- 最大那条 `20261011_122407`：5 段/5 簇，合并 96.15 → 求和 **0.08**。
  摊开看：5 段里 4 段单独测不出来（找不到两端静止 / 残差模型敏感 / 尾部静止不够 0.20 s），
  被跳过 ⟹ 求和只剩一段能测的。
- `20261011_122158`：2 段/1 簇，合并 3.92 → 求和 0.00（并簇后窗口变了）。
- `20261011_122416` / `122239`：合并与求和**都不行**，引擎自己就拒答了，
  这两条不构成逐簇的反例。
- **反过来看，`--fast` 那 6 条干净旧采集（1004xx×3 + 1437xx×2）逐簇求和几乎不变**
  （`|差|` 中位 0.00、最大 0.09 cm）—— 说明它只在"多段且有段测不出来"时才崩，
  不是普适地坏。这既是它还能救的理由，也是必须配判据的理由。

**判词**：逐簇求和把 "大声的错"（96 cm，明摆着离谱）变成了
"安静的错"（0.08 cm，看着体面）—— 正是本项目明确禁止的那类数。
**在有一个能把"敲击空段"与"被切碎的连续滑动"分开的判据之前，不许进引擎。**

用法：
    python sim/exp_split_risk.py            # 全库
    python sim/exp_split_risk.py --fast     # 只跑有真值的 6 条 + 122416/122239/122303
"""

import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2
from exp_real import load, ROOT, HARD, SOFT, TRUTH_CM
from exp_split import clusters_of, GAP_CLUSTER_S, rms

FOCUS = ["20261011_122416", "20261011_122239", "20261011_122303"]


def split_sum(acc, gyr, dn, up, dt, gap_min_s):
    m = measure2.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                         head_r=None, gate_motion=True)
    if not m["ok"]:
        return None
    merged = m["distanceM"] * 100
    total, nfail = 0.0, 0
    cls = clusters_of(m["groups"], dt, gap_min_s)
    for cl in cls:
        a, b = cl[0][0], cl[-1][1]
        mm = measure2.measure(acc, gyr, a, b, dt=dt, merge_mode="press",
                              head_r=None, gate_motion=True)
        if not mm["ok"]:
            nfail += 1
            continue
        total += mm["distanceM"] * 100
    return merged, total, m["motionS"], len(m["groups"]), len(cls), nfail


def one(nm, dt_cache=None):
    acc, gyr, dn, up, dt = load(os.path.join(ROOT, nm))
    return split_sum(acc, gyr, dn, up, dt, GAP_CLUSTER_S)


def main():
    fast = "--fast" in sys.argv
    names = sorted(HARD | SOFT) + FOCUS
    if not fast:
        names += [d for d in os.listdir(ROOT)
                  if os.path.isdir(os.path.join(ROOT, d))
                  and d not in set(names)]
        names = sorted(set(names))

    print("合并 vs 逐簇求和 —— 全库证伪扫描（GAP_CLUSTER_S = %.2f s）\n" % GAP_CLUSTER_S)
    print("%-18s %6s %4s %4s %4s | %9s %9s %9s %8s"
          % ("采集", "跨度s", "段", "簇", "拒", "合并cm", "求和cm", "差", "真值"))
    print("-" * 92)

    diffs, rows = [], []
    for nm in names:
        if not os.path.isdir(os.path.join(ROOT, nm)):
            continue
        try:
            r = one(nm)
        except Exception as e:
            print("%-18s  跳过（%s）" % (nm, e))
            continue
        if r is None:
            continue
        merged, total, span, ng, ncl, nf = r
        t = TRUTH_CM.get(nm)
        d = total - merged
        diffs.append(abs(d))
        rows.append((nm, span, ng, ncl, nf, merged, total, d, t))
        print("%-18s %6.2f %4d %4d %4d | %9.2f %9.2f %+9.2f %8s"
              % (nm, span, ng, ncl, nf, merged, total, d,
                 ("%.1f" % t) if t else "—"))

    print("-" * 92)
    ds = sorted(diffs)
    print("\n|求和−合并| 分布（n=%d）：中位 %.2f  均值 %.2f  P90 %.2f  最大 %.2f cm"
          % (len(ds), statistics.median(ds), statistics.fmean(ds),
             ds[int(len(ds) * 0.9)] if ds else 0, ds[-1] if ds else 0))

    big = sorted([r for r in rows if abs(r[7]) > 2.0], key=lambda r: -abs(r[7]))
    if big:
        print("\n|差| > 2 cm 的 %d 条（**逐条看是哪种工况，不下结论**）：" % len(big))
        for r in big[:20]:
            print("   %-18s 跨度%.2f  %d段/%d簇  合并 %7.2f → 求和 %7.2f  (%+.2f)"
                  % (r[0], r[1], r[2], r[3], r[5], r[6], r[7]))

    with_truth = [r for r in rows if r[8]]
    if with_truth:
        print("\n有真值的 %d 条（真值 15 cm）：" % len(with_truth))
        for r in with_truth:
            tag = "硬" if r[0] in HARD else "估"
            print("   %-18s(%s) 合并 %6.2f → 求和 %6.2f   （真值 %.1f）"
                  % (r[0], tag, r[5], r[6], r[8]))
        print("   合并 |误差| RMS %.2f / 求和 %.2f cm"
              % (rms([r[5] - r[8] for r in with_truth]),
                 rms([r[6] - r[8] for r in with_truth])))

    print("\n⚠ 143 系列三条的终点锚点本身就在动（比值 0.90~1.01），")
    print("  它们上面谁更准都说明不了什么 —— 列出来只是看逐簇会不会**变糟**。")


if __name__ == "__main__":
    main()
