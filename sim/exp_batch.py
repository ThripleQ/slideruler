"""实验 16：批量估计（v5）vs 双锚点去趋势（v4）—— 在**有精确真值**的仿真器上判对错。

为什么要先在仿真里做：真机只有 6 条真值、且真值本身 ±3%，分不清"算法变好"还是
"碰巧"。仿真里真值精确、噪声/游走可控，能直接回答"这个方法的收益是真的吗"。

判据：
  · 无噪声：v5 必须精确复原真值（验证实现正确，不是"看起来差不多"）；
  · 加噪声/零偏游走：v5 的误差应当**单调地**小于 v4，且随 T 的增长更缓；
  · 松手时仍在动：v5 不许出现 v4 那种"把真实位移扣掉"的塌陷。
"""
import argparse
import math
import statistics
import sys

sys.path.insert(0, '.')
from simgen import SimCfg, scenario, gen, down_i_of
import measure2 as M
import batch as B

CLEAN = dict(acc_noise=0.0, gyr_noise=0.0, acc_bias=(0, 0, 0), gyr_bias=(0, 0, 0),
             acc_scale=(0, 0, 0), acc_misalign=0.0)
# K40 量级：白噪声 + 常量零偏 + 零偏游走（游走才是主宰项）
NOISY = dict(acc_noise=0.010, gyr_noise=0.002, acc_bias=(0.020, -0.010, 0.030),
             gyr_bias=(0.002, -0.001, 0.003), acc_scale=(0.011, 0.011, 0.011),
             acc_misalign=0.004, gyr_bias_rw=0.02)
R_BODY = (0.0, -0.09, 0.0)


def one(name, kw, cfg_kw, seed=1, use_r=True, frame="body"):
    spec = scenario(name, cfg=SimCfg(**cfg_kw), r_body=R_BODY, **kw)
    spec["cfg"].seed = seed
    r = gen(spec["cfg"], spec["head_pos"], spec["att"], spec["r_body"], spec["t_end"])
    dt = 1 / spec["cfg"].hz
    dn = down_i_of(r["acc_ts"], spec["tA"] - 0.25)
    up = down_i_of(r["acc_ts"], spec["tB"] + 0.25)
    hr = R_BODY if use_r else None
    m4 = M.measure(r["acc"], r["gyr"], dn, up, dt=dt, merge_mode="press",
                   head_r=hr, gate_motion=False)
    m5 = B.measure_batch(r["acc"], r["gyr"], dn, up, dt=dt, head_r=hr, frame=frame)
    return spec["truth"], m4, m5


def cell(m4, m5, truth):
    def g(m):
        if m is None or not m.get("ok"):
            return None
        return m["distanceM"]
    return g(m4), g(m5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clean", action="store_true")
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--frame", default="body", choices=("body", "nav"))
    a = ap.parse_args()
    cfg = CLEAN if a.clean else NOISY
    tag = "无噪声（验证实现）" if a.clean else "K40 量级噪声 + 零偏游走 0.02 rad/s/√s"
    global FRAME
    FRAME = a.frame
    print("### %s   倾斜误差参数化坐标系: %s" % (tag, a.frame))

    cases = [
        ("slide_T0.5", "slide", dict(T=5.0, d=0.15, move_s=0.5), 15.0),
        ("slide_T0.7", "slide", dict(T=5.0, d=0.15, move_s=0.7), 15.0),
        ("slide_T1.0", "slide", dict(T=5.0, d=0.15, move_s=1.0), 15.0),
        ("slide_T1.4", "slide", dict(T=5.0, d=0.15, move_s=1.4), 15.0),
        ("slide_T2.0", "slide", dict(T=5.0, d=0.15, move_s=2.0), 15.0),
        ("slide_rot_T1.4", "slide_rotate", dict(T=5.0, d=0.15, move_s=1.4, rot_deg=15), 15.0),
        ("slide_wob_T1.4", "slide_wobble", dict(T=5.0, d=0.15, move_s=1.4, amp_deg=(3, 3, 4),
                                                drift_deg_per_s=(0, 0, 2)), 15.0),
        ("static", "static", dict(T=5.0), 0.0),
        ("static_rotate30", "static_rotate", dict(T=5.0, rot_deg=30), 0.0),
    ]
    print("%-20s %7s | %10s %10s | %10s %10s | %10s" % (
        "场景", "真值", "v4均值cm", "v5均值cm", "v4 σcm", "v5 σcm", "v5−v4中位"))
    print("-" * 92)
    for name, sname, kw, truth_cm in cases:
        v4s, v5s = [], []
        for s in range(a.n):
            truth, m4, m5 = one(sname, dict(kw), cfg, seed=s + 1, frame=FRAME)
            d4, d5 = cell(m4, m5, truth)
            if d4 is not None:
                v4s.append(d4 * 100)
            if d5 is not None:
                v5s.append(d5 * 100)
        f = lambda xs: ("%10.3f" % statistics.fmean(xs)) if xs else "%10s" % "-"
        sd = lambda xs: ("%10.3f" % statistics.pstdev(xs)) if len(xs) > 1 else "%10s" % "-"
        md = ("%+10.3f" % statistics.median([b - x for b, x in zip(v5s, v4s)])) if v4s and v5s else "-"
        print("%-20s %7.1f | %s %s | %s %s | %s" % (
            name, truth_cm, f(v4s), f(v5s), sd(v4s), sd(v5s), md))
    print("-" * 92)
    print("相对误差（有真值的场景）：")
    print("%-20s %7s | %9s %9s" % ("场景", "真值", "v4", "v5"))
    for name, sname, kw, truth_cm in cases:
        if truth_cm == 0:
            continue
        v4s, v5s = [], []
        for s in range(a.n):
            truth, m4, m5 = one(sname, dict(kw), cfg, seed=s + 1, frame=FRAME)
            d4, d5 = cell(m4, m5, truth)
            if d4 is not None:
                v4s.append((d4 * 100 - truth_cm) / truth_cm * 100)
            if d5 is not None:
                v5s.append((d5 * 100 - truth_cm) / truth_cm * 100)
        print("%-20s %7.1f | %8.1f%% %8.1f%%" % (
            name, truth_cm,
            statistics.fmean(v4s) if v4s else float("nan"),
            statistics.fmean(v5s) if v5s else float("nan")))


if __name__ == "__main__":
    main()
