"""实验：松手时手机还在动（早松手）—— 找停稳点能不能把它救回来。

真实起因：20261011_091811。用户按住 1.81 s 后松手，但那一刻手机还有 +1.2 cm/s 的
前向速度，松手后又往前走了 228 ms 才停住见顶（15.62 cm），随后被手拉回来。
按钮窗 [DOWN, UP] 因此少算了一段位移。

仿真里没法让"真值运动"在松手后继续，所以做法是：造一条完整的 15 cm 滑动，
然后把按钮窗人为截在运动中途（模拟早松手），看：
  A. 直接量（不处理）      → 误差随截断位置单调恶化，这就是用户能看见的失效
  B. 松手判据 + 顺延停稳点 → 应该把截掉的那段补回来，且与截断位置无关

用法: python exp_release.py
"""

import math
import statistics

import pipeline as P
from pipeline import G, Params, vdot, vnorm
import exp_slow as ES

DT = 0.0095138          # K40 实测
MARGIN_S = 0.3


def one(D, T, frac, seed, dt=DT):
    """把按钮窗截在运动进度的 frac 处，返回 (直接量, 顺延后, 真实终点位移)。"""
    cfg = ES.make_cfg(dt)
    tr = ES.gen_truth_2(cfg, D, T, seed=seed)
    n, n_pre, n_mot = tr["N"], tr["n_pre"], tr["n_mot"]
    prof = P.build_profile(tr)
    p = Params(static_mode="adaptive", assumed_dt=dt)
    qs, _ = P.solve_attitude(tr, prof.gyro_bias, prof.gravity_ref, p, dt=dt)
    a_lin = P.remove_gravity(tr, qs, g_mag=prof.gravity_mag)

    m = int(MARGIN_S / dt)
    lo = max(0, n_pre - m)
    hi_btn = min(n, n_pre + max(2, int(frac * n_mot)))
    hi_true = n_pre + n_mot

    cal = range(max(0, lo - int(2.0 / dt)), lo)
    b = math.sqrt(sum(statistics.fmean(a_lin[i][x] for i in cal) ** 2 for x in range(3)))

    def dist(hi):
        d, v, used, thr = P.integrate_cond_zupt3d(
            a_lin, lo, hi, dt, bias_bound=b, sigma=prof.acc_sigma, k_sigma=5.0)
        return d, abs(v), thr

    naive, v_end, thr = dist(hi_btn)
    end = hi_btn
    if v_end > thr:                     # 与条件 ZUPT 共用的松手判据
        end = P.stop_endpoint(a_lin, lo, hi_btn, dt=dt, thr=thr,
                              max_ext_s=1.2, dwell_s=0.05)
    fixed, _, _ = dist(end)
    truth_at_end, _, _ = dist(hi_true)
    return naive, fixed, truth_at_end, end - hi_btn


def main():
    print("=" * 96)
    print("早松手实验：15 cm 滑动，把按钮窗截在运动中途（模拟松手时手机还在动）")
    print(f"dt = {DT*1000:.4f} ms（K40 实测）   20 个种子   窗口左端 = DOWN - {MARGIN_S*1000:.0f} ms")
    print("=" * 96)
    print(f"{'截在运动':>9}{'':>2}{'直接量 中位/最大误差':>24}{'顺延后 中位/最大误差':>24}"
          f"{'顺延触发率':>12}{'平均顺延ms':>12}")
    for frac in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        naives, fixeds, exts = [], [], []
        for s in range(20):
            naive, fixed, truth_at_end, ext = one(0.15, 1.81, frac, s)
            naives.append(abs(naive / 0.15 - 1) * 100)
            fixeds.append(abs(fixed / 0.15 - 1) * 100)
            exts.append(ext)
        naives.sort()
        fixeds.sort()
        trig = sum(1 for e in exts if e > 0) / len(exts) * 100
        print(f"{frac*100:8.0f}%  "
              f"{statistics.median(naives):10.2f}%/{naives[-1]:6.2f}%      "
              f"{statistics.median(fixeds):10.2f}%/{fixeds[-1]:6.2f}%      "
              f"{trig:9.0f}%   {statistics.fmean(e*DT*1000 for e in exts):10.0f}")

    print()
    print("说明：")
    print("  · 直接量一栏的错误随截断位置单调恶化 —— 这就是按钮边界'没有平台'的样子，")
    print("    也就是用户能当场看见、但算法没法事后补救的失效。")
    print("  · 顺延一栏把终点推到'停稳区里前向位移最大'处，中位误差 ~1%，")
    print("    且与截断位置基本无关 —— 这一条把早松手从'整条作废'降成'正常误差'。")
    print("  · 顺延触发率 = 松手判据正确识别出'还在动'的比例（截断时应该 100%）。")
    print("  · frac=1.0 是'松手时刚好走完'，本来就不该顺延：这一行验证不误触发。")


if __name__ == "__main__":
    main()
