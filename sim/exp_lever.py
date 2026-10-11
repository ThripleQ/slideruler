"""杠杆臂 r 能不能从现有数据里辨识出来？

r 是"IMU 到测量头"的机体系偏移。它**不能**从设备资料查到（厂商不给），
也不该拍脑袋定 —— 拍错方向会让修正反向，把对的改成错的。

这里做一个直接的可辨识性实验：
  三条真值相同（都是用户估的 15 cm）但姿态变化各不相同的采集，
  如果用对了 r，它们修正后应当**收敛到同一个数**；用错了就会互相散开。
  于是"三条结果的离散度"就是 r 的目标函数 —— 扫 r 找它最小的位置。

注意这依赖真值相等（不要求真值精确），所以 15 cm 那 ±0.5 cm 反而无关紧要。
但 (R1−R0) 在三条里必须指向不同方向，否则某些分量根本不可辨识 —— 脚本会一并报告。
"""
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pipeline import qrot, vadd, vnorm, vsub
import measure2
from exp_real import load, ROOT

# 真值相同的采集（用户肉眼估的 15 cm ×3）
GROUP = ["20261011_100444", "20261011_100459", "20261011_100529"]


def collect():
    """跑一次算法，拿到每条的 IMU 位移矢量与两端姿态。"""
    out = []
    for nm in GROUP:
        d = load(os.path.join(ROOT, nm))
        if d is None:
            continue
        acc, gyr, dn, up, dt = d
        m = measure2.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                             head_r=None, gate_motion=True)
        if not m["ok"]:
            print("%s 拒答：%s" % (nm, m["reason"]))
            continue
        out.append((nm, m["pImuVec"], m["q0"], m["q1"]))
    return out


def head_dist(item, r):
    _, pI, q0, q1 = item
    return vnorm(vadd(pI, vsub(qrot(q1, r), qrot(q0, r))))


def main():
    items = collect()
    if len(items) < 2:
        print("可用采集不足，无法辨识")
        return

    # ---- 先看 (R1−R0) 的方向矩阵是否够"张开" ----
    print("每条的姿态差（机体系三根轴在导航系里的位移，单位 cm/m）：")
    basis = []
    for nm, pI, q0, q1 in items:
        ex = vsub(qrot(q1, (1, 0, 0)), qrot(q0, (1, 0, 0)))
        ey = vsub(qrot(q1, (0, 1, 0)), qrot(q0, (0, 1, 0)))
        ez = vsub(qrot(q1, (0, 0, 1)), qrot(q0, (0, 0, 1)))
        basis.append((nm, ex, ey, ez))
        print("  %s  ΔR·x̂=%s  ΔR·ŷ=%s  ΔR·ẑ=%s"
              % (nm, tuple(round(c, 4) for c in ex),
                 tuple(round(c, 4) for c in ey), tuple(round(c, 4) for c in ez)))

    # ---- 扫 r ----
    print("\n扫描 r = (0, -L, 0)：三条修正后的一致性")
    print("%8s %10s %10s %10s %10s %10s" % ("L(cm)", "444", "459", "529", "均值cm", "极差cm"))
    best = None
    for L in [x / 1000.0 for x in range(0, 141, 10)]:
        r = (0.0, -L, 0.0)
        ds = [head_dist(it, r) * 100 for it in items]
        spread = max(ds) - min(ds)
        print("%8.1f %10.2f %10.2f %10.2f %10.3f %10.3f"
              % (L * 100, ds[0], ds[1], ds[2], statistics.fmean(ds), spread))
        if best is None or spread < best[1]:
            best = (L, spread, statistics.fmean(ds))
    print("\n最一致：L=%.1f cm  极差 %.3f cm  均值 %.2f cm（用户估的 15 cm）"
          % (best[0] * 100, best[1] * 100, best[2]))

    # ---- 三轴各自扫一遍，看哪个方向有信息 ----
    print("\n各轴单独扫描（其余为 0）的极差最小值：")
    for axis, name in ((0, "x"), (1, "y"), (2, "z")):
        bb = None
        for L in [x / 1000.0 for x in range(-120, 121, 5)]:
            r = [0.0, 0.0, 0.0]
            r[axis] = L / 100.0
            ds = [head_dist(it, tuple(r)) * 100 for it in items]
            sp = max(ds) - min(ds)
            if bb is None or sp < bb[1]:
                bb = (L, sp)
        print("  %s 轴：最佳 %+.1f cm  极差 %.3f cm" % (name, bb[0], bb[1]))


if __name__ == "__main__":
    main()
