"""真机 41 条回归 —— 任何改动都必须过这一关。

仿真能证明"数学推导和实现对不对"，证明不了"真机的噪声是不是这个模型"
（ChatGPT 也提醒过这点）。所以仿真上再漂亮，回真机退化就是退化。

真值：分三档，别混（见下方 TRUTH_CM 的注释）。
最有用的那档是 2026-10-16 手操 20 cm 的 4 条：真值本身 ±1 cm（粗），
但**终点锚点终于干净**（点按式 + 滑完按住不动）。单列脚本 `sim/exp_20cm.py`。
其余没有真值的采集只能看**自洽性**：平台宽度、模型敏感性、拒答率有没有变差。

用法：
    python sim/exp_real.py                 # v4 全开
    python sim/exp_real.py --v3            # v3 基线
    python sim/exp_real.py --lever 0.09    # 指定杠杆臂长度
"""
import argparse
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                    "android", "captures", "slideruler")

# 真值（cm）。三档，**别混为一谈** —— 前两档是"尺子准不准"，第三档是"数据干不干净"：
#   HARD —— 2026-10-11 14:37 用 15 cm 直尺测的三次，用户明确"姿态比较随意"。
#           ⚠ 但这三条的**终点锚点在动**（|V_anchor|/v_peak = 0.90~1.01，见
#             sim/exp_anchor2.py）：标为"静止"的那 0.143 s 里积分速度等于运动峰值
#             速度，"松手"那一刻手机还在滑。⟹ **三条全部不可用**，不是精密基准。
#             保留在册只为让历史回归的数字不被改写。
#   SOFT —— 用户肉眼估的 15 cm，本身 ±0.5 cm（±3.3%），只能看趋势。
#   HAND20 —— 2026-10-16 23:05~23:07，换 20 cm 尺子重采、**整场手操**。
#             用户自报操作精度 ±1 cm（±5%）—— 真值本身比上面两档都粗。
#             但它换来了这一年里唯一"数据可信"的一批：终点锚点干净
#             （松手后静止 0.29~3.88 s，锚点漂移 0.012~0.072 m/s），
#             起点锚点也干净（角噪声 0.0002~0.0006 rad/s，桌面受控档）。
#             有效 4 条 = 用户划的 A 组里排除 230616（操作失误，引擎也独立拒答）。
#             前 4 条（230336/230414/230434/230452）同组但用户判"没练好"，不入册。
HARD = {"20261011_143719", "20261011_143725", "20261011_143734"}
SOFT = {"20261011_100444", "20261011_100459", "20261011_100529"}
HAND20 = {"20261016_230538", "20261016_230558",
          "20261016_230634", "20261016_230651"}
TRUTH_CM = {nm: 15.0 for nm in HARD | SOFT}
TRUTH_CM.update({nm: 20.0 for nm in HAND20})


def tier_of(nm):
    """真值分档标签。有三档，永远不要把它们平均到一起。"""
    if nm in HARD:
        return "硬"
    if nm in HAND20:
        return "手操20"
    return "估"


def load(d):
    acc_ts, acc = measure2.read_pairs(os.path.join(d, "acc.csv"))
    gyr_ts, gyr = measure2.read_pairs(os.path.join(d, "gyr.csv"))
    meta = open(os.path.join(d, "meta.txt"), encoding="utf-8").read()
    m = re.search(r"DOWN\s+event=\d+ ms\s+recv_uptime=\d+ ms\s+recv_elapsed=(\d+) ns", meta)
    if not m:
        return None
    dn = int(m.group(1))
    up = int(re.search(r"UP\s+event=\d+ ms\s+recv_uptime=\d+ ms\s+recv_elapsed=(\d+) ns",
                       meta).group(1))
    n = min(len(acc), len(gyr))
    acc, gyr, acc_ts = acc[:n], gyr[:n], acc_ts[:n]
    dt = statistics.median([acc_ts[i + 1] - acc_ts[i] for i in range(n - 1)]) / 1e9

    def locate(t):
        lo, hi = 0, n - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if acc_ts[mid] < t:
                lo = mid + 1
            else:
                hi = mid
        return lo
    return acc, gyr, locate(dn), locate(up), dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--v3", action="store_true", help="两项都关，回到 v3")
    ap.add_argument("--nolever", action="store_true")
    ap.add_argument("--nogate", action="store_true")
    ap.add_argument("--lever", type=float, default=0.09)
    ap.add_argument("--mode", default="press")
    a = ap.parse_args()
    lever = not (a.v3 or a.nolever)
    gate = not (a.v3 or a.nogate)
    r_vec = (0.0, -a.lever, 0.0)

    print("杠杆臂 %s  运动段禁修正 %s  模式 %s"
          % ("开(%.0fcm)" % (a.lever * 100) if lever else "关", "开" if gate else "关", a.mode))
    print("%-16s %9s %9s %8s %6s %5s  %s"
          % ("采集", "测量头cm", "IMUcm", "杠杆臂cm", "平台cm", "敏感", "备注"))
    print("-" * 104)

    names = sorted(x for x in os.listdir(ROOT) if os.path.isdir(os.path.join(ROOT, x)))
    rows = []
    n_ok = n_rej = 0
    for nm in names:
        d = load(os.path.join(ROOT, nm))
        if d is None:
            continue
        acc, gyr, dn, up, dt = d
        m = measure2.measure(acc, gyr, dn, up, dt=dt, merge_mode=a.mode,
                             head_r=r_vec if lever else None, gate_motion=gate)
        if not m["ok"]:
            n_rej += 1
            print("%-16s %9s %9s %8s %6s %5s  %s"
                  % (nm, "—", "—", "—", "—", "—", "拒：" + (m["reason"] or "")[:44]))
            continue
        n_ok += 1
        v = m["distanceM"] * 100
        t = TRUTH_CM.get(nm)
        tier = tier_of(nm)
        note = ""
        if t:
            note = "真值%.1f(%s) %+.1f%%" % (t, tier, (v - t) / t * 100)
        if m.get("modelSensitive"):
            note += " ⚠模型敏感"
        rows.append((nm, v, m.get("imuDistM", 0) * 100, m.get("leverM", 0) * 100,
                     m["plateau"], m.get("modelSensitive", False), t))
        print("%-16s %9.2f %9.2f %8.2f %6.1f~%.1f %5s  %s"
              % (nm, v, m.get("imuDistM", 0) * 100, m.get("leverM", 0) * 100,
                 m["plateau"][0], m["plateau"][1],
                 "是" if m.get("modelSensitive") else "", note))

    print("-" * 104)
    print("给出数值 %d / 拒答 %d" % (n_ok, n_rej))
    for tag, tier, absunc in (("硬真值(锚点在动,不可用)", "硬", 0.5),
                              ("肉眼估", "估", 0.5),
                              ("手操20cm(±1cm)", "手操20", 1.0)):
        sel = [r for r in rows if r[6] and tier_of(r[0]) == tier]
        if not sel:
            continue
        errs = [(r[1] - r[6]) for r in sel]
        rel = [e / r[6] * 100 for e, r in zip(errs, sel)]
        print("%s %d 条：相对误差 %s  平均 %+.2f%%  绝对误差中位 %+.2f cm（真值自身 ±%.1f cm）"
              % (tag, len(rel), "  ".join("%+.1f%%" % e for e in rel),
                 statistics.fmean(rel), statistics.median(errs), absunc))
    if rows:
        print("平台宽度中位数 %.2f cm" % statistics.median(
            [r[4][1] - r[4][0] for r in rows]))
        print("模型敏感 %d 条" % sum(1 for r in rows if r[5]))


if __name__ == "__main__":
    main()
