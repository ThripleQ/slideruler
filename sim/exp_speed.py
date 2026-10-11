"""实验 19：滑动时长 → 精度。（T³ 定律的定量标定）

硬真值三条（15 cm）的误差预算显示，主宰项是**陀螺零偏不稳定性**经姿态漂移
变成的重力泄漏：

    b_g ~ 0.012~0.020 rad/s（分段均值之差估计，静止段内都抹不掉）
    姿态漂移 θ = b_g·T      泄漏 = g·sin θ ≈ g·b_g·T
    位移误差 = ½·泄漏·T² = ½·g·b_g·T³      ← **∝ T³**

实测预算：T=0.69 s → 3.2 cm；T=1.44 s → 18 cm；T=1.49 s → 25 cm。
这是器件的物理极限，算法改不动 —— 但 T 是立方级杠杆：T 减半，误差降到 1/8。

本实验用仿真器（有精确真值，不靠肉眼估）标定这条曲线，回答：
  「该引导用户滑多快？」

用法：
    python sim/exp_speed.py            # 默认 K40 量级噪声，扫 0.25~2.0 s
    python sim/exp_speed.py --clean    # 无噪声，验证实现正确性
"""
import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pipeline import vnorm
import measure2 as M
from simgen import SimCfg, scenario, gen, down_i_of

# K40 / lsm6dso 量级（与 sim/exp_real.py 的 NOISY 一致）
NOISY = dict(acc_noise=0.010, gyr_noise=0.002,
             acc_bias=(0.020, -0.010, 0.030), gyr_bias=(0.002, -0.001, 0.003),
             acc_scale=(0.011, 0.011, 0.011), acc_misalign=0.004)
CLEAN = dict(acc_noise=0.0, gyr_noise=0.0, acc_bias=(0, 0, 0), gyr_bias=(0, 0, 0),
             acc_scale=(0, 0, 0), acc_misalign=0.0)


def one(T, d=0.15, cfg=None, seed=1, r_body=(0.0, 0.0, 0.0), gyr_rw=None):
    # move_s 就是滑动持续时间（scenario 里 tA..tB）
    c = dict(cfg or NOISY)
    if gyr_rw is not None:
        c["gyr_bias_rw"] = gyr_rw
    spec = scenario("slide", cfg=SimCfg(**c), r_body=r_body,
                    T=5.0, d=d, move_s=T)
    spec["cfg"].seed = seed
    r = gen(spec["cfg"], spec["head_pos"], spec["att"], spec["r_body"], spec["t_end"])
    dt = 1 / spec["cfg"].hz
    dn = down_i_of(r["acc_ts"], spec["tA"] - 0.25)
    up = down_i_of(r["acc_ts"], spec["tB"] + 0.25)
    m = M.measure(r["acc"], r["gyr"], dn, up, dt=dt, merge_mode="press",
                  head_r=r_body)
    return (m["distanceM"] if m["ok"] else None), spec["truth"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clean", action="store_true")
    ap.add_argument("--n", type=int, default=25)
    ap.add_argument("--d", type=float, default=0.15)
    ap.add_argument("--gyr-rw", type=float, default=None,
                    help="陀螺零偏随机游走 rad/s/√s（真机实测 ~0.02~0.03）")
    a = ap.parse_args()
    cfg = CLEAN if a.clean else NOISY

    Ts = [0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0]
    print("滑动距离 %.0f cm，%s，每组 %d 次，陀螺零偏游走 %s"
          % (a.d * 100, "无噪声" if a.clean else "K40 量级噪声", a.n,
             a.gyr_rw if a.gyr_rw is not None else "0（未建模！）"))
    print("%7s %7s %9s %9s %9s %9s %8s" %
          ("滑动T s", "真值cm", "均值cm", "偏差cm", "σ cm", "P90误差", "相对%"))
    print("-" * 68)
    curve = []
    for T in Ts:
        vals = []
        for s in range(a.n):
            v, truth = one(T, a.d, cfg, seed=s + 1, gyr_rw=a.gyr_rw)
            if v is not None:
                vals.append(v * 100)
        if not vals:
            print("%7.2f   全部拒答" % T)
            continue
        t = truth * 100
        errs = [v - t for v in vals]
        ae = sorted(abs(e) for e in errs)
        p90 = ae[int(0.9 * (len(ae) - 1))]
        curve.append((T, statistics.fmean(errs), p90))
        print("%7.2f %7.1f %9.2f %9.2f %9.2f %9.2f %8.1f" %
              (T, t, statistics.fmean(vals), statistics.fmean(errs),
               statistics.pstdev(vals), p90, abs(statistics.fmean(errs)) / t * 100))

    print("-" * 68)
    if len(curve) >= 2:
        print("T³ 检验（误差 / T³ 应当近似常数）：")
        for T, e, p90 in curve:
            print("   T=%.2f s  T³=%.3f  平均偏差 %+.2f cm  比值 %+.2f"
                  % (T, T ** 3, e, e / T ** 3 if T else 0))
        print()
        print("结论：把 T 从 1.5 s 压到 0.5 s，T³ 缩小 %.0f 倍"
              % (1.5 ** 3 / 0.5 ** 3))


if __name__ == "__main__":
    main()
