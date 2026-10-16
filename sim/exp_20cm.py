#!/usr/bin/env python3
"""
2026-10-16 夜 手操 20 cm 批次 —— 本项目唯一「数据本身可信 **且** 有尺子真值」的采集。
共 26 条（A 组 4 条 + 夜批 22 条）。用户 2026-10-17 确认：**真值全部 20 cm、
全部平放桌面推**。

为什么它值钱（前三批真值各自的毛病）：

  批 次                 真值来源            真值精度   数据本身
  --------------------  ------------------  ---------  --------------------------
  20261011_1004xx ×3    用户肉眼估 15 cm     ±0.5 cm   终点锚点干净，但真值粗
  20261011_1437xx ×3    15 cm 直尺 ×3        ±0.5 cm   **终点锚点在动**（0.90~1.01），废
  20261016 夜 ×26       20 cm 尺子、手操     **±1 cm**  **终点锚点终于干净**

配方 = 点按式边界 + 滑完手指继续按住不动 1~4 s。终点锚点漂移从 143 系列的
0.072~0.40 m/s 掉到 0.008~0.090 m/s，起点锚点角噪声 0.0002~0.0006 rad/s（桌面受控档）。

⚠ 真值自身的 ±1 cm（±5%）与引擎 RMS 同量级。所以这一批能回答的是
  「引擎有没有系统跑偏」，回答不了「引擎有没有 2 mm 级的偏差」。
  任何比 1 cm 更细的结论，都不能从这 26 条推出来。

⚠ 结论**只适用于平放桌面推**：这 26 条全部如此。手持悬空推时起止姿态变化
  大得多（见 exp_rot.py），外推会错。

用法：
    python sim/exp_20cm.py
"""
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2
from exp_real import HAND20, ROOT, TRUTH_CM, load

# 与 exp_real 默认一致（杠杆臂开、运动段禁修正）。这 4 条里 leverM 只有 0.01~0.04 cm，
# 开关它不改变结论 —— 滑动本身的转角太小，杠杆臂项够不着。
HEAD_R = (0.0, -0.09, 0.0)


def _pearson(xs, ys):
    """Pearson 相关系数。样本太少时不给"相关"——这里只在 n>=5 时调用。"""
    n = len(xs)
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0.0 or sy == 0.0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)


def main():
    truths = {TRUTH_CM[nm] for nm in HAND20}
    assert len(truths) == 1, "HAND20 的真值必须统一，否则不能整体汇总"
    truth = truths.pop()
    print("真值 %.0f cm（用户手操、自报 ±1 cm、全部平放桌面）  有效 %d 条\n"
          % (truth, len(HAND20)))
    print("%-16s %7s %4s %9s %9s %8s %8s"
          % ("采集", "运动s", "段", "裸积分cm", "引擎cm", "裸误差", "引擎误差"))
    print("-" * 72)

    rows = []                      # (nm, motionS, 裸误差, 引擎误差)
    for nm in sorted(HAND20):
        d = load(os.path.join(ROOT, nm))
        acc, gyr, dn, up, dt = d
        m = measure2.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                             head_r=HEAD_R, gate_motion=True)
        if not m["ok"]:
            print("%-16s %7s %4s 拒答：%s" % (nm, "—", "—", m["reason"]))
            continue
        bare = m["bareM"] * 100
        v = m["distanceM"] * 100
        be, ee = bare - truth, v - truth
        rows.append((nm, m["motionS"], be, ee))
        print("%-16s %7.2f %4d %9.2f %9.2f %+8.2f %+8.2f"
              % (nm, m["motionS"], m["nMerged"], bare, v, be, ee))

    print("-" * 72)
    if not rows:
        return
    bare_err = [r[2] for r in rows]
    eng_err = [r[3] for r in rows]
    n = len(rows)
    rms = lambda xs: (sum(x * x for x in xs) / len(xs)) ** 0.5
    mae = lambda xs: statistics.fmean(abs(x) for x in xs)

    print("裸积分  误差 RMS %.2f cm   平均 %+.2f cm   最大 |%.2f| cm"
          % (rms(bare_err), statistics.fmean(bare_err),
             max(abs(x) for x in bare_err)))
    print("引擎    误差 RMS %.2f cm   平均 %+.2f cm   最大 |%.2f| cm"
          % (rms(eng_err), statistics.fmean(eng_err),
             max(abs(x) for x in eng_err)))

    n_bare_pos = sum(1 for x in bare_err if x > 0)
    n_improve = sum(1 for be, ee in zip(bare_err, eng_err) if abs(ee) < abs(be))
    worse = [r for r in rows if abs(r[3]) >= abs(r[2])]
    print("\n裸积分误差为正 %d/%d 条 —— **不是恒为正**：22 条夜批里漂移有正有负，"
          % (n_bare_pos, n))
    print("  A 组那 4 条「全为正」是小样本巧合，不是规律。")
    print("引擎改善 %d/%d 条；平均 |误差| %.2f → %.2f cm；RMS 改良 %.2f×"
          % (n_improve, n, mae(bare_err), mae(eng_err),
             rms(bare_err) / rms(eng_err)))
    print("⚠ 整体 RMS 只是一个数，**掩盖单条成败**。引擎把本来准的拉偏的有 %d 条："
          % len(worse))
    for nm, _, be, ee in sorted(worse, key=lambda r: -abs(r[3]))[:5]:
        print("     %s  裸 %+.2f → 引擎 %+.2f" % (nm, be, ee))
    if worse:
        print("  被拉偏条目运动时长中位 %.2f s vs 全体 %.2f s"
              % (statistics.median(r[1] for r in worse),
                 statistics.median(r[1] for r in rows)))
    print("  引擎不是单条的保险 —— 整体均值在修，**单条可能更差**。")

    r1 = _pearson([r[1] for r in rows], [abs(x) for x in eng_err])
    print("\n窗口-误差：|引擎误差| ~ 运动时长 Pearson %.3f"
          % (r1))
    print("带符号引擎误差 ~ 运动时长 Pearson %.3f"
          % _pearson([r[1] for r in rows], eng_err))
    print("按运动时长三分位看带符号误差中位：")
    by_t = sorted(rows, key=lambda r: r[1])
    k = max(1, n // 3)
    for lo, hi, tag in ((0, k, "短"), (k, 2 * k, "中"), (2 * k, n, "长")):
        seg = by_t[lo:hi]
        if seg:
            print("  %s %.2f~%.2f s  n=%2d  中位 %+6.2f cm"
                  % (tag, seg[0][1], seg[-1][1], len(seg),
                     statistics.median(r[3] for r in seg)))
    print("\n⚠ 真值自身 ±1.0 cm 与引擎 RMS %.2f cm 同量级 ⟹ 只能说「均值上没跑偏」"
          % rms(eng_err))
    print("  （平均 %+.2f cm），**不能说「误差 %.2f cm」**——那会把真值的精度当成引擎的。"
          % (statistics.fmean(eng_err), rms(eng_err)))


if __name__ == "__main__":
    main()
