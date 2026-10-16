"""2026-10-17 凌晨这批（8 条）为什么比昨晚 26 条脏那么多 —— 只用数据说话。

用户自述「操作比较准确、真值 20 cm」，但引擎读数散布在 10~49 cm，3 条拒答。
本脚本**不用真值**，只比较两批的「数据自身健康度」指标，看退化是操作带进来的
还是算法退化了。

判据（全部来自当次数据，不依赖真值）：
  · 终点锚点比 = anchorDrift / vPeak  —— 松手后那段"静止"里还剩多少速度
  · 起点角噪声 startTremorRadS        —— 起手那段手机放没放稳
  · σ_a（加速度计噪声，meta 里）
  · 运动段结构（段数 / 在动时长 / 空白占比）
  · 修正量占比 = modelGap / distance   —— 答案有多少来自残差模型

用法：
    python sim/exp_20261017.py
"""

import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2
from exp_real import load, ROOT, HAND20, TRUTH_CM

NEW = [
    "20261017_004039", "20261017_004058", "20261017_004107", "20261017_004111",
    "20261017_004117", "20261017_004135", "20261017_004148", "20261017_004156",
]


def mode_of(meta):
    m = re.search(r"滑动测量（(点按式|按住式)）", meta)
    return m.group(1) if m else "—"


def sigma_a(meta):
    m = re.search(r"σ_a\s+([\d.]+)", meta)
    return float(m.group(1)) if m else None


def per_seg(acc, gyr, dt, groups):
    """逐段单独测（按钮窗 = 该段，零余量），看每段自己带多少位移。"""
    out = []
    for a, b in groups:
        mm = measure2.measure(acc, gyr, a, b, dt=dt, merge_mode="press",
                              head_r=None, gate_motion=True)
        out.append(((b - a) * dt,
                    mm["distanceM"] * 100 if mm["ok"] else None,
                    (mm.get("reason") or "")[:20]))
    return out


def main():
    print("2026-10-17 凌晨批（8 条）逐条体检 —— 真值若均为 20 cm\n")
    for nm in NEW:
        d = os.path.join(ROOT, nm)
        meta = open(os.path.join(d, "meta.txt"), encoding="utf-8").read()
        ld = load(d)
        if ld is None:
            print("%-16s %s  （无按钮边界）" % (nm, mode_of(meta)))
            continue
        acc, gyr, dn, up, dt = ld
        m = measure2.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                             head_r=None, gate_motion=True)
        print("=" * 96)
        print("%s   %s   σ_a=%s   边界跨度 %.2f s"
              % (nm, mode_of(meta), sigma_a(meta), (up - dn) * dt))
        if not m["ok"]:
            print("  拒答：%s" % m["reason"])
            # 拒答也看段结构（groups 仍给了）
            if m.get("groups"):
                for dur, dist, why in per_seg(acc, gyr, dt, m["groups"]):
                    print("    段 %.2f s → %s" % (dur, ("%.2f cm" % dist)
                                                  if dist is not None else "拒:" + why))
            continue
        drift = m.get("anchorDrift", 0.0)
        vpk = m.get("vPeakMps", 1e-9)
        print("  读数 %6.2f cm   裸积分 %7.2f   修正量 %6.2f cm（答案的 %.0f%%）"
              % (m["distanceM"] * 100, m.get("bareM", 0) * 100,
                 m.get("modelGapM", 0) * 100,
                 m.get("modelGapM", 0) / max(m["distanceM"], 1e-9) * 100))
        print("  在动 %.2f s / 窗口 %.2f s   段数 %d   终点静止 %.2f s"
              % (m["motionS"], (m["window"][1] - m["window"][0]) * dt,
                 len(m["groups"]), m.get("tailStaticS", 0.0)))
        print("  终点锚点 %5.3f m/s  峰值 %5.3f m/s   比值 %.2f %s"
              % (drift, vpk, drift / max(vpk, 1e-9),
                 "⚠ 锚点在动" if drift / max(vpk, 1e-9) > 0.5 else ""))
        print("  起点角噪声 %.4f rad/s（τ=%.2f s）   平台 %.2f~%.2f cm"
              % (m.get("startTremorRadS", 0.0), m.get("startTremorTauS", 0.0),
                 m["plateau"][0], m["plateau"][1]))
        for dur, dist, why in per_seg(acc, gyr, dt, m["groups"]):
            print("    段 %.2f s → %s" % (dur, ("%.2f cm" % dist)
                                          if dist is not None else "拒:" + why))

    # ---- 两批对比（只用数据，不用真值）----
    print("\n" + "=" * 96)
    print("两批「数据健康度」对比（全部指标不依赖真值）\n")

    def health(names, tag):
        rows = []
        for nm in sorted(names):
            ld = load(os.path.join(ROOT, nm))
            if ld is None:
                continue
            acc, gyr, dn, up, dt = ld
            m = measure2.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                                 head_r=None, gate_motion=True)
            if not m["ok"]:
                continue
            meta = open(os.path.join(ROOT, nm, "meta.txt"), encoding="utf-8").read()
            sa = sigma_a(meta)
            rows.append({
                "nm": nm,
                "mode": mode_of(meta),
                "ratio": m.get("anchorDrift", 0) / max(m.get("vPeakMps", 1e-9), 1e-9),
                "tremor": m.get("startTremorRadS", 0.0),
                "sigma": sa if sa is not None else float("nan"),
                "gap": m.get("modelGapM", 0) / max(m["distanceM"], 1e-9),
                "blank": 1.0 - m["motionS"] / max(
                    (m["window"][1] - m["window"][0]) * dt, 1e-9),
                "nseg": len(m["groups"]),
            })
        if not rows:
            return
        med = lambda k: statistics.median([r[k] for r in rows
                                           if r[k] == r[k]])
        print("%s  %d 条出数" % (tag, len(rows)))
        print("  终点锚点比 中位 %.2f   最大 %.2f   >0.5 的 %d 条"
              % (med("ratio"), max(r["ratio"] for r in rows),
                 sum(1 for r in rows if r["ratio"] > 0.5)))
        print("  起点角噪声 中位 %.4f   最大 %.4f"
              % (med("tremor"), max(r["tremor"] for r in rows)))
        print("  σ_a        中位 %.3f" % med("sigma"))
        print("  修正量占比 中位 %.0f%%" % (med("gap") * 100))
        print("  窗口空白比 中位 %.0f%%   段数中位 %.0f"
              % (med("blank") * 100, med("nseg")))
        modes = {}
        for r in rows:
            modes[r["mode"]] = modes.get(r["mode"], 0) + 1
        print("  按钮模式：%s\n" % "  ".join("%s×%d" % (k, v)
                                            for k, v in sorted(modes.items())))

    health(HAND20, "【昨晚 26 条】")
    health([n for n in NEW], "【今晨 8 条】")


if __name__ == "__main__":
    main()
