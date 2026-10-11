"""定位实验：积分区间与锚点选择如何决定精度（判断误差的真正来源）。

真值：滑动段 = 采样点 [200, 350]，位移恰好 0.8 m，v(200) = v(350) = 0。
"""

import math
import statistics

import pipeline as P
from pipeline import DT, D_TRUE, Params, SensorCfg, gen_truth

IDEAL = dict(acc_bias=0, acc_noise=0, gyro_bias=0, gyro_noise=0, gyro_instab=0,
             wobble_deg=0, tilt_deg=0)
NOMINAL = dict()


def integrate(a1, lo, hi):
    v = [0.0]
    p = [0.0]
    for i in range(lo + 1, hi + 1):
        v.append(v[-1] + (a1[i - 1] + a1[i]) / 2 * DT)
        p.append(p[-1] + (v[-2] + v[-1]) / 2 * DT)
    return v, p


def trial(a1, lo, hi, drift=True):
    v, p = integrate(a1, lo, hi)
    if not drift:
        return abs(p[-1])
    da = v[-1] / ((hi - lo) * DT)
    a1c = list(a1)
    for i in range(lo + 1, hi + 1):
        a1c[i] = a1[i] - da
    return abs(integrate(a1c, lo, hi)[1][-1])


def build(cfg_kw, seed=7):
    cfg = SensorCfg(**cfg_kw)
    truth = gen_truth(cfg, seed=seed)
    qs, _ = P.solve_attitude(truth, *P.calibrate(truth, truth["n_pre"])[:2], Params())
    a_lin = P.remove_gravity(truth, qs)
    u, _, _ = P.pca_axis(a_lin[201:349])
    return truth, [sum(a_lin[i][k] * u[k] for k in range(3)) for i in range(len(a_lin))]


def show(label, d):
    print(f"  {label:<40} {d * 100:7.2f} cm   误差 {(d - D_TRUE) / D_TRUE * 100:+7.2f}%")


def block(title, kw):
    print("=" * 100)
    print(title)
    print("=" * 100)
    truth, a1 = build(kw)
    show("真值窗口 [200,350]", trial(a1, 200, 350))
    show("检出窗口 [200,378]（滑窗检测延迟 0.28s）", trial(a1, 200, 378))
    show("记录起点→检出终点 [0,378]", trial(a1, 0, 378))
    print()

    print(f"  {'起点偏移(终点固定350)':<26}{'距离':>10}{'误差':>10}")
    for k in [0, 1, 2, 5, 10, 30]:
        d = trial(a1, 200 + k, 350)
        print(f"  晚 {k:>2} 采样 = {k * DT * 1000:>4.0f} ms{'':<6}{d * 100:>10.2f}"
              f"{(d - D_TRUE) / D_TRUE * 100:>9.2f}%")
    print()
    print(f"  {'终点偏移(起点固定200)':<26}{'距离':>10}{'误差':>10}")
    for k in [0, 1, 2, 5, 10, 30]:
        d = trial(a1, 200, 350 + k)
        print(f"  晚 {k:>2} 采样 = {k * DT * 1000:>4.0f} ms{'':<6}{d * 100:>10.2f}"
              f"{(d - D_TRUE) / D_TRUE * 100:>9.2f}%")
    print()


def main():
    block("实验组 1：理想传感器 + 真值姿态 + 理想主轴（隔离一切传感器误差）", IDEAL)
    block("实验组 2：标称误差全开（零偏/噪声/手部转动，姿态用滤波器估计）", NOMINAL)

    print("=" * 100)
    print("实验组 3：起点晚了半个窗长时，各误差源各自的贡献（终点固定 350）")
    print("=" * 100)
    for name, kw in [("理想传感器", IDEAL), ("标称传感器", NOMINAL)]:
        truth, a1 = build(kw)
        show(f"{name}，起点 200", trial(a1, 200, 350))
        show(f"{name}，起点 215（晚 150 ms）", trial(a1, 215, 350))


if __name__ == "__main__":
    main()
