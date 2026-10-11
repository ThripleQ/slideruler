"""杠杆臂标定：能不能用"测量头顶住不动、手机绕它转"这个动作反解出 r？

背景
----
exp_lever.py 证明了 r **辨识不出来**：三条真值相同的真机采集，
姿态差的方向太单一，r 的贡献和噪声/其他误差混在一起（极差最小点就在 r=0）。
所以不能靠"多测几次拟合出来"，必须设计一个**信噪比足够高**的专用动作。

原理
----
把测量头顶在一个固定点上不动，只让手机绕它转。
测量头的位移在物理上恒为 0，于是刚体关系直接变成对 r 的线性方程：

    Δp_H = Δp_I + [R(t1) − R(t0)]·r = 0
    ⇒   [R(t1) − R(t0)]·r = −Δp_I

左边是 3×3 矩阵（可以从姿态估出来），右边是 3 向量（可以从 IMU 积分出来）。
一次转动给 3 个方程，绕不同轴转几次就能最小二乘解出 r 的 3 个分量。

信噪比是这个方法成立的关键：转 90° 时 |ΔR·r| = 2|r|·sin(45°) ≈ 1.41|r|，
r=6 cm 就是 8.5 cm 的信号；而 ZUPT 之后 Δp_I 的误差在 1 cm 量级 —— 够了。
对比 ordinary 滑动：转 15° 只有 0.26|r| ≈ 1.6 cm 的信号，被噪声埋掉，所以辨识不了。

本脚本在仿真里验证：解出来的 r 能否恢复真值、误差多大、需要几次转动。
"""
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pipeline import qrot, qmul, qnorm, q_from_euler, vadd, vsub, vnorm
from simgen import SimCfg, gen, down_i_of, _slerp, att_const, smoothstep5
import measure2


def solve_ls3(A, b, lam=1e-12):
    """最小二乘解 A x = b（A 是若干 3×3 矩阵，b 是若干 3 向量）。正规方程 AᵀA x = Aᵀb。

    注意 ΔR = R(t1) − R(t0) **永远奇异**：任何旋转矩阵都有特征值 1（旋转轴），
    所以 (RΔ − I) 在旋转轴方向上是零空间。一次转动只能定出 r 垂直于该轴的分量，
    沿轴的分量根本不可观测 —— 必须绕**不同的轴**多转几次才能把三个分量都解出来。
    这正是 exp_lever.py 里"三条普通滑动辨识不出 r"的同一件事的另一种表现。
    """
    M = [[0.0] * 4 for _ in range(3)]
    for Ai, bi in zip(A, b):
        for p in range(3):
            for q in range(3):
                M[p][q] += (Ai[0][p] * Ai[0][q] + Ai[1][p] * Ai[1][q]
                            + Ai[2][p] * Ai[2][q])
            M[p][3] += Ai[0][p] * bi[0] + Ai[1][p] * bi[1] + Ai[2][p] * bi[2]
    for i in range(3):
        M[i][i] += lam
    for c in range(3):
        piv = max(range(c, 3), key=lambda r_: abs(M[r_][c]))
        if abs(M[piv][c]) < 1e-14:
            return None
        M[c], M[piv] = M[piv], M[c]
        for r_ in range(3):
            if r_ == c:
                continue
            f = M[r_][c] / M[c][c]
            for k in range(c, 4):
                M[r_][k] -= f * M[c][k]
    return tuple(M[i][3] / M[i][i] for i in range(3))


def pivot_rotation(r_true, axis_euler, deg, cfg, seed):
    """一次"支点转动"：测量头始终在原点，手机绕它转 deg 度。"""
    q0 = qnorm((1.0, 0.0, 0.0, 0.0))
    dq = q_from_euler(*[math.radians(deg) * c for c in axis_euler])
    q1 = qnorm(qmul(q0, dq))
    tA, tB, tEnd = 1.0, 2.0, 3.2
    # 必须用 smoothstep 而不是线性插值：线性插值让**角速度在起止点突变**
    # （0 → 常数 → 0），二阶差分在跳变处给出 v/h ≈ 47 m/s² 的假冲击，
    # 积分出来比真值大 4 倍。真实的转动是平滑加速再减速的。
    u = lambda t: smoothstep5((t - tA) / (tB - tA))
    att = lambda t: _slerp(q0, q1, u(t))
    head = lambda t: (0.0, 0.0, 0.0)      # 测量头钉死不动
    cfg.seed = seed
    r = gen(cfg, head, att, r_true, tEnd)
    dt = 1.0 / cfg.hz
    dn = down_i_of(r["acc_ts"], tA - 0.25)
    up = down_i_of(r["acc_ts"], tB + 0.25)
    m = measure2.measure(r["acc"], r["gyr"], dn, up, dt=dt, merge_mode="press",
                         head_r=None, gate_motion=True)
    if not m["ok"]:
        return None
    # [R(t1) − R(t0)] 的三列
    A = [tuple(vsub(qrot(m["q1"], e), qrot(m["q0"], e)))
         for e in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))]
    A = [[A[0][i], A[1][i], A[2][i]] for i in range(3)]   # 转成行=方程的形式
    b = tuple(-c for c in m["pImuVec"])
    return A, b, m


# 绕不同轴的转动（机体系欧拉角方向）
ROTATIONS = [
    ((0, 0, 1), 90, "偏航 90°"),
    ((1, 0, 0), 90, "横滚 90°"),
    ((0, 1, 0), 90, "俯仰 90°"),
    ((1, 1, 0), 70, "斜轴 70°"),
    ((0, 0, 1), -90, "偏航 -90°"),
    ((1, 0, 1), 80, "斜轴 80°"),
]


def main():
    r_true = (0.0, -0.060, 0.0)     # 假设真值：测量头在 IMU 的 -y 方向 6 cm
    print("真值 r = %s m" % (tuple(round(c, 4) for c in r_true),))

    for label, cfg in [("无噪声", SimCfg(acc_noise=0.0, gyr_noise=0.0,
                                         acc_bias=(0, 0, 0), gyr_bias=(0, 0, 0))),
                       ("K40 量级噪声", SimCfg(acc_noise=0.010, gyr_noise=0.002,
                                               acc_bias=(0.020, -0.010, 0.030),
                                               gyr_bias=(0.002, -0.001, 0.003)))]:
        A_all, b_all = [], []
        print("\n=== %s ===" % label)
        print("%-12s %10s %10s  %s" % ("转动", "信号cm", "残差cm", "解出的 r (cm)"))
        for k, (axis, deg, name) in enumerate(ROTATIONS):
            out = pivot_rotation(r_true, axis, deg, cfg, seed=k + 1)
            if out is None:
                print("%-12s  拒答" % name)
                continue
            A, b, m = out
            sig = vnorm(b) * 100
            A_all.append(A)
            b_all.append(b)
            r_cum = solve_ls3(A_all, b_all)
            if r_cum is None:
                print("%-12s %10.2f %10s  （秩不足，还需绕别的轴转）" % (name, sig, "—"))
                continue
            err = vnorm(vsub(r_cum, r_true)) * 100
            print("%-12s %10.2f %10.2f  (%+.2f, %+.2f, %+.2f)"
                  % (name, sig, err, *[r_cum[i] * 100 for i in range(3)]))
        if A_all:
            r_f = solve_ls3(A_all, b_all)
            print("最终 r = (%+.2f, %+.2f, %+.2f) cm，真值 (%+.2f, %+.2f, %+.2f) cm，"
                  "误差 %.2f cm"
                  % tuple(list(r_f[i] * 100 for i in range(3))
                          + list(r_true[i] * 100 for i in range(3))
                          + [vnorm(vsub(r_f, r_true)) * 100]))


if __name__ == "__main__":
    main()
