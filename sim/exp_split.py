"""窗长是被「点按敲击」撑大的 —— 逐簇测量再求和。

## 发现（2026-10-17，26 条手操 20 cm 平放桌面）

把每次采集的运动段**逐段单独测**（按钮窗收紧到该段，零余量）之后，结构一目了然：

    234859  DOWN=1.08  UP=6.28
      0.97-1.17 (0.20s) →  0.01 cm     ← 敲击
      3.79-4.81 (1.03s) → 20.11 cm     ← 真正的滑动
      6.16-6.42 (0.26s) →  0.01 cm     ← 敲击
      6.82-7.01 (0.19s) → 拒答          ← 敲击（尾部静止不够）
      合并（现役）27.81 cm   逐段求和 20.12 cm   真值 20

**每次"运动"里只有一段带着位移，其余全是 0.01~0.03 cm 的空段**，而且它们的时间位置
死死贴着 DOWN / UP 两个按钮时刻 —— 是**手指敲击屏幕的机械脉冲**被静止检测认成了运动。
26 条里几乎每条的首尾各有一到两个这样的敲击段。

## 为什么这会毁掉读数

现役 `merge_mode="press"` 取"按钮按住期间所有运动段的**并集**"，于是窗口从
敲击段一直跨到敲击段：

    234859  敲击到敲击 5.45 s      （真正在动 1.67 s）
    234423  敲击到敲击 4.79 s      （真正在动 2.16 s）
    230538  敲击到敲击 6.75 s      （真正在动 2.22 s）

而双锚点的修正量约等于 **T·(V₀+V₁)/2**，T 是**窗口跨度**（推导见下）。
停顿/敲击对真实位移零贡献，却把 T 撑大 2~3 倍 ⟹ 修正量同步放大 ⟹ §5.8 那条
"≤3.8 s 全准、>3.8 s 全崩"的分界线。**分界的物理来源不是"滑得快不快"，是"窗口里
混进了多少不属于这次滑动的空白时间"。**

### 修正量为什么 ∝ T

`zupt_two_anchor` 解 e(t) = A·t + B 使两侧锚点段的**速度均值**为 0。设两锚点的平均
时刻 t₀≈0、t₁≈T，则 A ≈ (V₁−V₀)/T、B ≈ V₀，被扣掉的位移是

    ∫₀ᵀ (A·t + B) dt = A·T²/2 + B·T ≈ T·(V₀ + V₁)/2

**与真实运动多久无关，只与窗口跨度成正比。** 缩短 T 是唯一能直接压下它的手段。

## 本脚本的做法

1. 按**间隙**把运动段聚类：间隙 < GAP_CLUSTER_S 的归为一簇 —— 它们之间没有足够静止，
   塞不下各自的锚点，物理上只能一起测；
2. 每簇**单独**走一次 `measure()`（按钮窗 = 该簇本身，零余量，保证只选中这一簇）；
3. 各簇距离**求和**。敲击段自己只有 0.01 cm，加进来无害，关键是**它们不再拉长窗口**。

### GAP_CLUSTER_S = 0.30 s 是先验的，不是挑出来的

窗口两侧各要 MARGIN_S = 0.15 s 的静止当锚点，一次间隙要同时喂饱左右两簇 ⟹
**2 × MARGIN_S = 0.30 s**。与"测量头/真值"无关。

（扫描曲线一并打印，但**不作为选参依据** —— 挑最优就是拿真值调参。）

## ⚠ 已证伪的失败模式（见 `exp_split_risk.py`，**这条阻止它进引擎**）

在旧语料上跑全库对照，抓到一条灾难性的：

    20261011_122407   5 段/5 簇   合并 96.15 cm → 求和 **0.08 cm**

它的 5 段里有 4 段**单独测不出来**（找不到两端静止 / 残差模型敏感 / 尾部静止不够），
被脚本跳过 ⟹ 求和只剩那一段能测的 0.08 cm。

**这是本项目最忌讳的那类错误：把"大声的错"变成"安静的错"。** 96 cm 明摆着离谱，
0.08 cm 看着体面 —— 而两者都错。（122416 / 122239 则是合并与求和都不行，
引擎自己就拒答了，那两条不构成反例。）

根因是 v3 那个旧伤：**静止检测分不出匀速和静止**，一次连续滑动会被切成好几段，
段间的"静止"其实是匀速运动。现役 `merge_mode="press"` 取并集正是为了挡它
（122416 一次按住被切 5 段，只取一段 → 0.59 cm）。逐簇求和把这个保护拆了。

⟹ **不能直接上。** 要上必须先有一个判据，把"敲击空段"和"被切碎的连续滑动"
区分开 —— 而这两者在**当次数据**上长得很像（都是 0.2~0.3 s 的短段，都不带位移）。

## 已知边界

- 求和用的是**标量**距离 Σ|Δp_i|。26 条都是沿尺子单向直推，方向一致，
  所以 Σ|Δp_i| ≈ |Σ Δp_i|；**有往返动作时会失效**，那时必须改成矢量求和。
- 簇太小（static run < STATIC_MIN_S）会被 `measure()` 拒答并**跳过** ⟹ 少算。
  本批 26 条里有 4 条（230558/234423/234735/234859）各有一个簇被拒，
  恰好都是敲击段，所以求和仍对 —— **但这是运气，不是保证**（122407 就是反例）。

## 结论强度

`[实测 N=26，in-sample；已在旧语料上证伪]`。判据不依赖真值（只看间隙与锚点可用性），
但"好不好"是用真值评的。**不要接进引擎** —— 见上面的失败模式。

用法：
    python sim/exp_split.py
"""

import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2
from exp_real import load, ROOT, HAND20, TRUTH_CM

MARGIN_S = 0.15           # 与 SlideMeasure.kt 一致：窗口 = 运动段 ± MARGIN_S
GAP_CLUSTER_S = 2 * MARGIN_S   # 先验：一次间隙要喂饱左右两簇的锚点余量


def clusters_of(groups, dt, gap_min_s):
    """把运动段按间隙聚类。间隙 = 上一段结束 → 下一段开始。"""
    if not groups:
        return []
    out = [[groups[0]]]
    for g in groups[1:]:
        if (g[0] - out[-1][-1][1]) * dt < gap_min_s:
            out[-1].append(g)
        else:
            out.append([g])
    return out


def measure_sum(acc, gyr, dn, up, dt, gap_min_s):
    """逐簇测量再求和。返回 (总距离 m, 簇数, 拒答数, 最长簇窗 s, 合并窗 s)。"""
    m = measure2.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                         head_r=None, gate_motion=True)
    if not m["ok"]:
        return None
    merged_s = (m["window"][1] - m["window"][0]) * dt
    total, nfail, worst = 0.0, 0, 0.0
    for cl in clusters_of(m["groups"], dt, gap_min_s):
        a, b = cl[0][0], cl[-1][1]
        # 按钮窗 = 该簇本身（零余量）：merge_mode="press" 只收与按钮窗重叠的段，
        # 这样能保证选中的就是这一簇。
        mm = measure2.measure(acc, gyr, a, b, dt=dt, merge_mode="press",
                              head_r=None, gate_motion=True)
        if not mm["ok"]:
            nfail += 1
            continue
        total += mm["distanceM"]
        worst = max(worst, (mm["window"][1] - mm["window"][0]) * dt)
    ncl = len(clusters_of(m["groups"], dt, gap_min_s))
    return total, ncl, nfail, worst, merged_s


def rms(xs):
    return (sum(x * x for x in xs) / len(xs)) ** 0.5


def main():
    do_scan = "--scan" in sys.argv
    print("逐簇测量求和 vs 合并长窗 —— 真值 20 cm，26 条手操（平放桌面）\n")
    print("GAP_CLUSTER_S = %.2f s（先验 = 2×MARGIN_S，不是从真值挑的）\n" % GAP_CLUSTER_S)

    rows = []
    for nm in sorted(HAND20):
        acc, gyr, dn, up, dt = load(os.path.join(ROOT, nm))
        m = measure2.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                             head_r=None, gate_motion=True)
        if not m["ok"]:
            continue
        r = measure_sum(acc, gyr, dn, up, dt, GAP_CLUSTER_S)
        if r is None:
            continue
        tot, ncl, nfail, worst, merged_s = r
        truth = TRUTH_CM[nm]
        busy = sum((b - a) * dt for a, b in m["groups"])
        rows.append((nm, m["motionS"], busy, len(m["groups"]), ncl, nfail,
                     merged_s, worst, m["distanceM"] * 100,
                     m["distanceM"] * 100 - truth, tot * 100, tot * 100 - truth))

    print("%-18s %6s %6s %3s %3s %3s | %6s %6s | %8s %8s | %8s %8s"
          % ("采集", "跨度s", "在动s", "段", "簇", "拒", "合并窗", "最长簇",
             "合并", "合并误", "求和", "求和误"))
    print("-" * 104)
    for r in rows:
        print("%-18s %6.2f %6.2f %3d %3d %3d | %6.2f %6.2f | %8.2f %+8.2f | %8.2f %+8.2f"
              % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11]))
    print("-" * 104)

    me = [r[9] for r in rows]
    se = [r[11] for r in rows]
    n = len(rows)
    for tag, es in (("合并长窗", me), ("逐簇求和", se)):
        print("%s  RMS %5.2f cm   平均 %+5.2f   最大 |%5.2f|   |误差|≤1cm %2d/%d"
              % (tag, rms(es), statistics.fmean(es), max(abs(x) for x in es),
                 sum(1 for x in es if abs(x) <= 1.0), n))
    b = sum(1 for r in rows if abs(r[11]) < abs(r[9]) - 1e-9)
    w = sum(1 for r in rows if abs(r[11]) > abs(r[9]) + 1e-9)
    print("逐簇更准 %d / 更差 %d / 持平 %d" % (b, w, n - b - w))
    print("\n窗口跨度：合并中位 %.2f s → 最长簇中位 %.2f s（真正的积分窗口）"
          % (statistics.median([r[6] for r in rows]),
             statistics.median([r[7] for r in rows])))

    tail = [r for r in rows if r[1] > 3.8]
    if tail:
        print("\n长窗那 %d 条（跨度 >3.8 s，合并方案的失败区）：" % len(tail))
        print("  合并 RMS %.2f → 逐簇 %.2f cm" % (rms([r[9] for r in tail]),
                                              rms([r[11] for r in tail])))
    short = [r for r in rows if r[1] <= 3.8]
    if short:
        print("短窗那 %d 条（合并本来就准，逐簇应当几乎不变）：" % len(short))
        print("  合并 RMS %.2f → 逐簇 %.2f cm" % (rms([r[9] for r in short]),
                                              rms([r[11] for r in short])))

    if not do_scan:
        print("\n（GAP 扫描要跑约 1000 次 measure()，默认跳过；想看加 `--scan`）")
        print("\n⚠ in-sample：26 条、真值 ±1 cm。判据不依赖真值（只看间隙），")
        print("  但「好不好」是真值评的。接引擎前必须：① 全新批次复核；")
        print("  ② Python/Kotlin 两侧同步 + 单测；③ 想清楚往返动作时怎么求和（见文件头边界）。")
        return

    print("\nGAP_CLUSTER_S 扫描（**只看形状，不用来挑参数** —— 挑最优=拿真值调参）：")
    print("  %10s %8s %8s %6s %6s" % ("GAP", "RMS", "平均", "簇数", "拒答"))
    for g in (0.16, 0.20, 0.30, 0.40, 0.60, 0.90, 1.50):
        es, ncl, nf = [], 0, 0
        for nm in sorted(HAND20):
            acc, gyr, dn, up, dt = load(os.path.join(ROOT, nm))
            r = measure_sum(acc, gyr, dn, up, dt, g)
            if r is None:
                continue
            ncl += r[1]
            nf += r[2]
            es.append(r[0] * 100 - TRUTH_CM[nm])
        print("  %10.2f %8.2f %+8.2f %6d %6d" % (g, rms(es), statistics.fmean(es), ncl, nf))

    print("\n⚠ in-sample：26 条、真值 ±1 cm。判据不依赖真值（只看间隙），")
    print("  但「好不好」是真值评的。接引擎前必须：① 全新批次复核；")
    print("  ② Python/Kotlin 两侧同步 + 单测；③ 想清楚往返动作时怎么求和（见文件头边界）。")


if __name__ == "__main__":
    main()
