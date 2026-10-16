"""实验 19：滑动时长 → 精度。

## ⚠ 本实验的立论已被证伪（2026-10-16 注，见 docs/v4-head-model.md §5.1）

文件名头原来写的是"主宰项是陀螺零偏不稳定性，位移误差 ∝ T³"。现在知道那不对：
库里 `20261011_085650`（76 s 纯静止录制）上，零偏 5 s 窗口均值全落在 ±3e-5 rad/s 内、
无单调漂移，Allan deviation 从 10 ms 到 9.7 s 一路单调下降、无平台也无回升
—— **这台 K40 的零偏不游走**。`gyro_bias_instability`/`attitude_drift_budget`
与界面上的"慢滑动"报警都已删除。

**本脚本保留但降级为参考**：
- 它跑的是**仿真器**，而仿真器自己的陀螺游走模型（`gyr_bias_rw=0.02`，无界游走）
  连形状都错了（真机单调下降，它单调上升）—— 重标之前它的结论不能当真（§6.3）；
- 它下面打印的"T³ 检验"里的 T³ 列是**已废弃的框架**，看"平均偏差"列就好；
- 想回答"滑多快更准"，正确做法是**重采硬真值**（§7 #1），不是在仿真里扫 T。

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
                    help="陀螺零偏随机游走 rad/s/√s。⚠ 曾经标的 0.02~0.03 是拿 P90 算的，"
                         "已被证伪 —— 真机不游走（76 s 静止录制，见文件头）")
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
