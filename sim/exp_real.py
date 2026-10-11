"""真机 41 条回归 —— 任何改动都必须过这一关。

仿真能证明"数学推导和实现对不对"，证明不了"真机的噪声是不是这个模型"
（ChatGPT 也提醒过这点）。所以仿真上再漂亮，回真机退化就是退化。

真值：只有 3 条是用户肉眼估的 15 cm（本身 ±0.5 cm），只能看趋势、不能当精密基准。
其余 38 条没有真值，只能看**自洽性**：平台宽度、模型敏感性、拒答率有没有变差。

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

# 真值（cm）。两档可信度，别混为一谈：
#   HARD —— 2026-10-11 14:37 用 15 cm 直尺测的三次，用户明确"姿态比较随意"。
#           这是本项目第一批**硬真值**，能判到 1% 量级。
#   SOFT —— 用户肉眼估的 15 cm，本身 ±0.5 cm（±3.3%），只能看趋势。
TRUTH_CM = {
    "20261011_143719": 15.0,
    "20261011_143725": 15.0,
    "20261011_143734": 15.0,
}
HARD = set(TRUTH_CM)
TRUTH_CM.update({
    "20261011_100444": 15.0,
    "20261011_100459": 15.0,
    "20261011_100529": 15.0,
})


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
        note = ""
        if t:
            note = ("真值%.1f(硬) %+.1f%%" if nm in HARD else "真值%.1f(估) %+.1f%%") % (
                t, (v - t) / t * 100)
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
    for tag, sel in (("硬真值", [r for r in rows if r[0] in HARD]),
                     ("肉眼估", [r for r in rows if r[6] and r[0] not in HARD])):
        if not sel:
            continue
        errs = [(r[1] - r[6]) / r[6] * 100 for r in sel]
        print("%s %d 条相对误差：%s  平均 %+.2f%%  最大 |%.2f%%|"
              % (tag, len(errs), "  ".join("%+.1f%%" % e for e in errs),
                 statistics.fmean(errs), max(abs(e) for e in errs)))
    if rows:
        print("平台宽度中位数 %.2f cm" % statistics.median(
            [r[4][1] - r[4][0] for r in rows]))
        print("模型敏感 %d 条" % sum(1 for r in rows if r[5]))


if __name__ == "__main__":
    main()
