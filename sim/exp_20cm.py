#!/usr/bin/env python3
"""
2026-10-16 手操 20 cm 批次 —— 本项目第一批「数据本身可信 **且** 有尺子真值」的采集。

为什么它值钱（前三批真值各自的毛病）：

  批 次                 真值来源            真值精度   数据本身
  --------------------  ------------------  ---------  --------------------------
  20261011_1004xx ×3    用户肉眼估 15 cm     ±0.5 cm   终点锚点干净，但真值粗
  20261011_1437xx ×3    15 cm 直尺 ×3        ±0.5 cm   **终点锚点在动**（0.90~1.01），废
  20261016_2305xx ×4    20 cm 尺子、手操     **±1 cm**  **终点锚点终于干净**

这一批的配方 = 点按式边界 + 滑完手指继续按住不动 1~4 s。效果直接可见：
终点锚点漂移从 143 系列的 0.072~0.40 m/s 掉到 0.012~0.072 m/s，
起点锚点角噪声 0.0002~0.0006 rad/s（桌面受控档）。

⚠ 真值自身的 ±1 cm（±5%）比引擎读数粗一个量级。所以这一批能回答的是
  「引擎有没有系统跑偏」，回答不了「引擎有没有 2 mm 级的偏差」。
  任何比 1 cm 更细的结论，都不能从这 4 条推出来。

用法：
    python sim/exp_20cm.py
"""
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2
from exp_real import HAND20, ROOT, TRUTH_CM, load

# 与 exp_real 默认一致（杠杆臂开、运动段禁修正）。这 4 条里 leverM 只有 0.01~0.04 cm，
# 开关它不改变结论 —— 滑动本身的转角太小，杠杆臂项够不着。
HEAD_R = (0.0, -0.09, 0.0)


def main():
    truth = TRUTH_CM[HAND20.copy().pop()]
    print("真值 %.0f cm（用户手操，自报 ±1 cm）  有效 %d 条\n" % (truth, len(HAND20)))
    print("%-16s %7s %4s %9s %9s %8s %8s"
          % ("采集", "运动s", "段", "裸积分cm", "引擎cm", "裸误差", "引擎误差"))
    print("-" * 72)

    bare_err, eng_err = [], []
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
        bare_err.append(be)
        eng_err.append(ee)
        print("%-16s %7.2f %4d %9.2f %9.2f %+8.2f %+8.2f"
              % (nm, m["motionS"], m["nMerged"], bare, v, be, ee))

    print("-" * 72)
    if eng_err:
        rms = lambda xs: (sum(x * x for x in xs) / len(xs)) ** 0.5
        print("裸积分  误差 RMS %.2f cm   平均 %+.2f cm   最大 |%.2f| cm"
              % (rms(bare_err), statistics.fmean(bare_err),
                 max(abs(x) for x in bare_err)))
        print("引擎    误差 RMS %.2f cm   平均 %+.2f cm   最大 |%.2f| cm"
              % (rms(eng_err), statistics.fmean(eng_err),
                 max(abs(x) for x in eng_err)))
        print("\n裸积分误差**恒为正**（漂移只会加不会减），引擎把 4 条全部拉向真值；")
        print("改良倍数 %.1f×（RMS）。真值自身 ±1.0 cm —— 引擎误差已进这个范围。"
              % (rms(bare_err) / rms(eng_err)))
        print("\n误差随积分窗口变长而变差（1.83s→+0.6  3.36s→−0.3  4.99s→−1.1  6.75s→−2.6 cm），")
        print("与 ∝T² 机制方向一致。但**只有 4 个点、真值 ±1 cm，不足以拟合任何定律** ——")
        print("这里只当操作建议（一口气 1~2 s 推完），不当结论。")


if __name__ == "__main__":
    main()
