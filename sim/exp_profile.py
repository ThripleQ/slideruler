"""决定性测试：起点起步平缓时，指南的滑窗静止检测会滞后多少、代价多大。"""

import math

import pipeline as P
from pipeline import DT, D_TRUE, Params, SensorCfg, gen_truth


def eval_window(a1, lo, hi):
    v = [0.0]
    p = [0.0]
    for i in range(lo + 1, hi + 1):
        v.append(v[-1] + (a1[i - 1] + a1[i]) / 2 * DT)
        p.append(p[-1] + (v[-2] + v[-1]) / 2 * DT)
    da = v[-1] / ((hi - lo) * DT)
    a1c = list(a1)
    for i in range(lo + 1, hi + 1):
        a1c[i] = a1[i] - da
    v2 = [0.0]
    for i in range(lo + 1, hi + 1):
        v2.append(v2[-1] + (a1c[i - 1] + a1c[i]) / 2 * DT)
    return abs(sum((v2[k - 1] + v2[k]) / 2 * DT for k in range(1, len(v2))))


def main():
    print("=" * 104)
    print("起步形态对静止检测的影响（0.8 m / 1.5 s / 100 Hz）")
    print("=" * 104)
    for profile, desc in [("sine", "正弦速度剖面：起步瞬间加速度即峰值（最友好）"),
                          ("cycloid", "摆线剖面：起步加速度为 0，缓慢起推（更像真实手）")]:
        cfg = SensorCfg(profile=profile)
        truth = gen_truth(cfg, seed=7)
        b_g, g_ref, _ = P.calibrate(truth, truth["n_pre"])
        qs, _ = P.solve_attitude(truth, b_g, g_ref, Params())
        a_lin = P.remove_gravity(truth, qs)
        u, _, _ = P.pca_axis(a_lin[201:349])
        a1 = [sum(a_lin[i][k] * u[k] for k in range(3)) for i in range(len(a_lin))]

        print(f"\n### {desc}")
        a_pk = max(abs(x) for x in truth["a_true"])
        print(f"  加速度峰值 {a_pk:.3f} m/s²"
              f"，起步后 0.05s 时 a = {truth['a_true'][205]:.3f} m/s²"
              f"，0.1s 时 {truth['a_true'][210]:.3f} m/s²")
        for mode, mdesc in [("guide_magnitude", "指南判据一 |‖f‖-9.81|<0.2"),
                            ("linear_accel", "‖a_lin‖<0.12"),
                            ("accel_variance", "窗内加速度 σ<0.05"),
                            ("oracle", "真值窗口")]:
            p = Params(static_mode=mode)
            flags = P.detect_static(truth, b_g, p, a_lin=a_lin)
            lo, hi = P.motion_window(flags, truth, p)
            d = eval_window(a1, lo, hi)
            lag = lo - 200
            print(f"  {mdesc:<28} 起点 {lo:>3}（滞后 {lag * 10:>3} ms） 终点 {hi:>3}  "
                  f"距离 {d * 100:6.2f} cm  误差 {(d - D_TRUE) / D_TRUE * 100:+7.2f}%")

    print()
    print("=" * 104)
    print("为什么起步平缓时要滞后这么久：判据一需要横向加速度到 ~1.99 m/s² 才越过阈值")
    print("=" * 104)
    G = 9.81
    a_need = math.sqrt((G + 0.2) ** 2 - G ** 2)
    print(f"  无零偏时，|‖f‖-9.81| 要达到 0.2 需要横向加速度 a = {a_need:.3f} m/s²")
    cfg = SensorCfg(profile="cycloid")
    t = gen_truth(cfg, seed=7)
    a_pk = max(abs(x) for x in t["a_true"])
    print(f"  摆线剖面峰值 {a_pk:.3f} m/s²，达到 {a_need:.3f} 需要 sin(2πτ/T) = "
          f"{a_need / a_pk:.3f} → τ = {math.asin(a_need / a_pk) / (2 * math.pi) * 1.5:.3f} s")
    print(f"  → 起步后 {math.asin(a_need / a_pk) / (2 * math.pi) * 1.5 * 1000:.0f} ms 才第一次触发，")


if __name__ == "__main__":
    main()
