"""v1 与 v2 在同一批真机采集上的逐条对照。

为什么要有这个脚本
------------------
v2 是"用户说按着不动也出距离"之后推倒重来的产物。要证明它真的治好了，光有仿真不够 ——
仿真里两端的静止段是纯静止，v1 的条件 ZUPT 反而表现不错（见 exp_v2.py 组 A）。
真实采集里那段静止带着**手引起的微姿态误差**，重力泄漏进水平通道，残余零偏远大于仿真，
v1 的自相矛盾判据就在那里显形。

本脚本把两个版本并排打在同一张表上，并把 v1 的内部依据（残余零偏、v_end、阈值、
它自己的审计结论）一并带出来 —— 因为 v1 的失效模式不是"报错"，而是
**审计全绿、却给出了 18.20 cm**。

v1 的脚本（`slide_pipeline.py`）保持冻结不改，所以这里用抓 stdout 的方式取值。
真值一栏只有用户目测报过的几次（15 cm ×3），其余留空。
"""
import io
import os
import contextlib
import re
import sys

import measure2 as M2
import slide_pipeline as SP


ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "android", "captures", "slideruler")

# 用户肉眼估过的真值（cm）。只有这几条能算误差。
TRUTH_CM = {
    "20261011_100444": 15.0,
    "20261011_100459": 15.0,
    "20261011_100529": 15.0,
}

RE_DIST = re.compile(r"★ 距离 = ([\d.]+) cm")
RE_RAW = re.compile(r"（三维裸积分 ([\d.]+) cm）")
RE_PATH = re.compile(r"→ (裸积分|去趋势)")
RE_VEND = re.compile(r"前向终点速度 v_end = ([+-][\d.]+) m/s\s+阈值 = ([\d.]+) m/s")
RE_BIAS = re.compile(r"残余零偏 \|b\| = ([\d.]+) m/s²")
RE_AUDIT = re.compile(r"(DOWN|UP)\s+±0\.3s\s+max‖a_lin‖ =\s+([\d.]+)\s+([✓✗])")
RE_CAL = re.compile(r"\|g\|=([\d.]+)\s+k̂=([\d.]+)\s+σ_a=([\d.]+)")


def run_v1(d):
    """跑 v1 管线并把它打印的东西抓回来（v1 只打印，没有返回值）。"""
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            SP.analyse(d)
    except Exception as e:  # v1 对退化输入没有防线
        return {"err": f"{type(e).__name__}: {e}"}
    s = buf.getvalue()
    out = {}
    if m := RE_DIST.search(s):
        out["d"] = float(m.group(1))
    if m := RE_RAW.search(s):
        out["raw"] = float(m.group(1))
    if m := RE_PATH.search(s):
        out["path"] = m.group(1)
    if m := RE_VEND.search(s):
        out["vend"] = float(m.group(1))
        out["thr"] = float(m.group(2))
    if m := RE_BIAS.search(s):
        out["bias"] = float(m.group(1))
    if m := RE_CAL.search(s):
        out["g"] = float(m.group(1))
        out["kcal"] = float(m.group(2))
        out["sa"] = float(m.group(3))
    aud = {k: (float(v), o) for k, v, o in RE_AUDIT.findall(s)}
    out["audit"] = aud
    return out


def fmt_v1(r):
    if "err" in r:
        return f"抛异常（{r['err'][:34]}）"
    if "d" not in r:
        return "没给距离"
    path = r.get("path", "?")
    return f"{r['d']:6.2f} cm  [{path}]"


def main():
    names = sorted(x for x in os.listdir(ROOT) if os.path.isdir(os.path.join(ROOT, x)))
    print("=" * 168)
    print("v1（条件 ZUPT）vs v2（两段静止夹运动）—— 同一批真机采集，逐条对照")
    print("-" * 168)
    print(f"{'采集':>16} {'真值':>6} | {'v1':^28} | {'v2':^30} | "
          f"{'v1 内部依据':^38}")
    print("-" * 168)

    for nm in names:
        d = os.path.join(ROOT, nm)
        truth = TRUTH_CM.get(nm)
        t_str = f"{truth:.0f}" if truth else "-"

        v1 = run_v1(d)
        v2 = M2.run_dir(d, verbose=False)

        # --- v1 这一列 ---
        c1 = fmt_v1(v1)
        if truth and "d" in v1:
            c1 += f" {((v1['d'] - truth) / truth * 100):+6.1f}%"

        # --- v2 这一列 ---
        if v2 is None:
            c2 = "无触摸标记"
        elif not v2["ok"]:
            c2 = f"拒答：{(v2.get('reason') or '')[:22]}"
        elif v2.get("noMotion"):
            c2 = "判没动 → 0.00 cm"
        else:
            c2 = f"{v2['distanceM']*100:6.2f} cm  平台 {v2['plateau'][0]:.2f}~{v2['plateau'][1]:.2f}"
            if truth:
                c2 += f" {((v2['distanceM']*100 - truth) / truth * 100):+5.1f}%"
            if v2.get("warn") or v2["tailStaticS"] < M2.MARGIN_S:
                c2 += " ⚠"

        # --- v1 的内部依据 ---
        if "err" in v1:
            c3 = "-"
        else:
            bits = []
            if "bias" in v1:
                bits.append(f"|b|={v1['bias']:.4f}")
            if "vend" in v1:
                bits.append(f"v_end={v1['vend']:+.4f}/{v1['thr']:.4f}")
            for k in ("DOWN", "UP"):
                if k in v1.get("audit", {}):
                    bits.append(f"{k[0]}{v1['audit'][k][1]}")
            c3 = " ".join(bits)

        print(f"{nm:>16} {t_str:>6} | {c1:<28} | {c2:<30} | {c3:<38}")

    print("-" * 168)
    print("说明：")
    print("  · v1 的 |b| 是「静止段残余零偏」，v_end/阈值 是条件 ZUPT 的判据 ——")
    print("    纯静止窗里 v_end 恒 ≈ |b|·T ≤ 阈值（阈值第一项就是 |b|·T），所以判据**必然**成立。")
    print("  · D/U 是 v1 自己的 ±0.3 s 巡检结论（✓ 静止）—— v1 在静止窗上巡检全绿，")
    print("    却给出十几到上百厘米；它的问题不是报错，是**自信地报错**。")
    print("  · v2 的 ⚠ 表示尾部静止不够（终点锚点没得用），这个数不可靠。")
    print("  · 真值只有用户目测的几条（±0.5 cm 量级），其余留空。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
