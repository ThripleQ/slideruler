"""诊断：误差到底由谁主导？逐项消融 + 残差测量 + 检测器对比。"""

import math
import statistics

import pipeline as P
from pipeline import (DT, D_TRUE, Params, SensorCfg, calibrate, gen_truth,
                      remove_gravity, run_pipeline, solve_attitude)


def ablate(label, cfg_kw=None, mode="oracle", perfect_attitude=False):
    cfg = SensorCfg(**(cfg_kw or {}))
    p = Params(static_mode=mode)
    truth = gen_truth(cfg, seed=7)
    if perfect_attitude:
        r = run_with_q(truth, p, truth["q_true"])
    else:
        r = run_pipeline(truth, cfg, p, seed=7)
    if not r.get("ok", True):
        print(f"  {label:<30} 失败：{r.get('reason')}")
        return r
    e = (r["dist"] - D_TRUE) / D_TRUE * 100
    print(f"  {label:<30} 窗口[{r['lo']:>3},{r['hi']:>3}]  距离 {r['dist'] * 100:7.2f} cm"
          f"  误差 {e:+7.2f}%")
    return r


def run_with_q(truth, p, qs):
    a_lin = remove_gravity(truth, qs)
    flags = P.detect_static(truth, (0.0, 0.0, 0.0), p, a_lin=a_lin)
    lo, hi = P.motion_window(flags, truth, p)
    u, _, _ = P.pca_axis(a_lin[lo:hi + 1])
    axis = u if p.use_pca else (1.0, 0.0, 0.0)
    a1 = [sum(a_lin[i][k] * axis[k] for k in range(3)) for i in range(len(a_lin))]
    v, pp = P.integrate(a1, lo, hi)
    da = v[-1] / ((hi - lo) * DT)
    a1c = [a1[i] - da if lo <= i <= hi else a1[i] for i in range(len(a1))]
    v2, p2 = P.integrate(a1c, lo, hi)
    return {"ok": True, "dist": abs(p2[-1]), "dist_raw": abs(pp[-1]),
            "lo": lo, "hi": hi, "v_err": v[-1]}


def main():
    print("=" * 100)
    print("诊断 1：误差源消融（窗口用真值，隔离检测器误差；只剩算法误差）")
    print("=" * 100)
    print("  真值滑动窗口应为 [201, 349]，长度 1.49 s\n")
    ablate("理想传感器（误差全关）", dict(acc_bias=0, acc_noise=0, gyro_bias=0,
                                     gyro_noise=0, gyro_instab=0, wobble_deg=0))
    ablate("只留加速度计零偏 0.05 m/s²", dict(acc_noise=0, gyro_bias=0, gyro_noise=0,
                                        gyro_instab=0, wobble_deg=0))
    ablate("只留加速度计噪声 0.02 m/s²", dict(acc_bias=0, gyro_bias=0, gyro_noise=0,
                                        gyro_instab=0, wobble_deg=0))
    ablate("只留陀螺噪声 0.002 rad/s", dict(acc_bias=0, acc_noise=0, gyro_bias=0,
                                       gyro_instab=0, wobble_deg=0))
    ablate("只留陀螺零偏 0.5 °/s", dict(acc_bias=0, acc_noise=0, gyro_noise=0,
                                     gyro_instab=0, wobble_deg=0))
    ablate("只留手部转动 1°", dict(acc_bias=0, acc_noise=0, gyro_bias=0,
                                gyro_noise=0, gyro_instab=0))
    ablate("全部打开（标称）")

    print()
    print("  对照：姿态换成真值（隔离互补滤波）")
    ablate("理想传感器 + 完美姿态", dict(acc_bias=0, acc_noise=0, gyro_bias=0,
                                    gyro_noise=0, gyro_instab=0, wobble_deg=0),
           perfect_attitude=True)
    ablate("标称误差 + 完美姿态", perfect_attitude=True)

    print()
    print("=" * 100)
    print("诊断 2：静止检测判据对比（同一段录制的滑动窗口识别结果）")
    print("=" * 100)
    print(f"  真值窗口 [201, 349]\n")
    cfg, p0 = SensorCfg(), Params()
    truth = gen_truth(cfg, seed=7)
    b_g, g_ref, _ = calibrate(truth, truth["n_pre"])
    qs, _ = solve_attitude(truth, b_g, g_ref, p0)
    a_lin = remove_gravity(truth, qs)
    for mode, desc in [("guide_magnitude", "指南原版 |‖f‖-9.81|<0.2"),
                       ("linear_accel", "去重力后 ‖a_lin‖<0.12"),
                       ("accel_variance", "窗内加速度标准差<0.05"),
                       ("oracle", "真值窗口（上界）")]:
        p = Params(static_mode=mode)
        flags = P.detect_static(truth, b_g, p, a_lin=a_lin)
        lo, hi = P.motion_window(flags, truth, p)
        seg = "无" if lo is None else f"[{lo}, {hi}] 长 {(hi - lo) * DT:.2f}s"
        print(f"  {desc:<30} → 检出窗口 {seg}")

    print()
    print("=" * 100)
    print("诊断 3：去重力后的残差（标称工况，滑动段真值 [201,349]）")
    print("=" * 100)
    print(f"  {'区间':<10}{'ax 残差均值':>14}{'ay 残差均值':>14}{'az 残差均值':>14}"
          f"{'ax 残差std':>13}")
    for name, lo, hi in [("静止前段", 0, 190), ("滑动段", 205, 345), ("静止后段", 400, truth["N"])]:
        ex = [a_lin[i][0] - truth["a_true"][i] for i in range(lo, hi)]
        print(f"  {name:<10}{statistics.fmean(ex):>14.5f}"
              f"{statistics.fmean(a_lin[i][1] for i in range(lo, hi)):>14.5f}"
              f"{statistics.fmean(a_lin[i][2] for i in range(lo, hi)):>14.5f}"
              f"{statistics.pstdev(ex):>13.5f}")
    a_perf = remove_gravity(truth, truth["q_true"])
    ex = [a_perf[i][0] - truth["a_true"][i] for i in range(205, 345)]
    print(f"  {'真值姿态':<10}{statistics.fmean(ex):>14.5f}"
          f"{statistics.fmean(a_perf[i][1] for i in range(205, 345)):>14.5f}"
          f"{statistics.fmean(a_perf[i][2] for i in range(205, 345)):>14.5f}"
          f"{statistics.pstdev(ex):>13.5f}")

    print()
    print("=" * 100)
    print("诊断 4：姿态估计质量与 α 的带宽")
    print("=" * 100)
    for lo, hi, name in [(0, 190, "静止前段"), (205, 345, "滑动段"), (400, truth["N"], "静止后段")]:
        errs = []
        for i in range(lo, hi):
            d = P.qnorm(P.qmul(P.qconj(truth["q_true"][i]), qs[i]))
            errs.append(math.degrees(2 * math.acos(min(1.0, abs(d[0])))))
        print(f"  {name:<8} 姿态角误差 均值 {statistics.fmean(errs):.4f}°  "
              f"最大 {max(errs):.4f}°")
    print(f"\n  α = 0.05 → 互补滤波时间常数 1/α = {1 / 0.05:.0f} s，"
          f"而记录总长仅 {truth['N'] * DT:.1f} s")


if __name__ == "__main__":
    main()
