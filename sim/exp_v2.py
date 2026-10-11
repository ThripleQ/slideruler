"""仿真对照：v2 链路（静止段夹窗口 + 双锚点 ZUPT + Mahony 参考）vs 现役链路。

真机只有两个肉眼估的「15 cm」标签，而且记录尾部常常没有静止段，判不出算法好坏。
这里用 K40 实测参数造带精确真值的记录，逐一注入失效模式：

  · 纯静止（D=0）        → 看会不会凭空造出距离（用户实际报的毛病）
  · 运动中缓慢旋转        → 看姿态参考系的敏感性（用户的「对姿态要求很严格」）
  · 按钮窗比运动长        → 看「按住不放」的静止段会不会被算进距离
  · 污染校准段            → 看校准失败时会不会给个数

用法: python exp_v2.py
"""
import math
import random
import statistics

import pipeline as P
from pipeline import G, Params, SensorCfg, vdot, vnorm, vscale
import exp_slow as ES
import measure2 as M2

DT = 0.0095138


# ------------------------------------------------------------------ v1 参考实现

def v1_distance(tr, down_i, up_i, dt, g_mag_mode="measured", cond=True):
    """现役链路：写死的校准窗（DOWN 前 2 s）+ PCA 一维投影 + 条件 ZUPT。"""
    n = tr["N"]
    prof = P.build_profile(tr)
    p = Params(static_mode="adaptive", assumed_dt=dt)
    qs, _ = P.solve_attitude(tr, prof.gyro_bias, prof.gravity_ref, p, dt=dt)
    a_lin = P.remove_gravity(tr, qs, g_mag=prof.gravity_mag if g_mag_mode == "measured" else G)
    lo, hi = down_i, up_i
    u, l1, l2 = P.pca_axis(a_lin[lo:hi + 1])
    axis = u if vdot(u, (1.0, 0.0, 0.0)) >= 0 else vscale(u, -1.0)
    a1 = [vdot(a_lin[i], axis) for i in range(n + 1)]
    cal = range(max(0, lo - int(2.0 / dt)), lo)
    sigma = statistics.pstdev([a1[i] for i in cal])
    b = abs(statistics.fmean(a1[i] for i in cal))
    if cond:
        d, _, _, _ = P.integrate_cond_zupt(a1, lo, hi, dt, bias_bound=b, sigma=sigma)
    else:
        _, pr = P.integrate(a1, lo, hi, dt=dt)
        d = abs(pr[-1])
    return d


def v1_attitude_off(tr, down_i, up_i, dt):
    """现役姿态（body 参考），其余同 v2 —— 用来单独隔离姿态参考系那一项。"""
    n = tr["N"]
    prof = P.build_profile(tr)
    p = Params(static_mode="adaptive", assumed_dt=dt)
    qs, _ = P.solve_attitude(tr, prof.gyro_bias, prof.gravity_ref, p, dt=dt)
    a_lin = P.remove_gravity(tr, qs, g_mag=prof.gravity_mag)
    # 用真值窗口 + 无条件去趋势，隔离出姿态那一项
    import measure2 as M
    T = (up_i - down_i) * dt
    v, pth = M.integrate_from(a_lin, down_i, up_i, dt)
    dr = [v[-1][x] / T for x in range(3)] if T > 1e-9 else [0, 0, 0]
    acc2 = [tuple(a_lin[i][x] - dr[x] for x in range(3)) for i in range(down_i, up_i + 1)]
    _, p2 = M.integrate_from(acc2, 0, up_i - down_i, dt)
    return vnorm(p2[-1])


# ------------------------------------------------------------------ 造记录

def make(D, T_mot, wobble_deg, seed, dt=DT, gyro_bias=0.0, acc_bias=0.092,
         press_early=0.25, release_late=0.30):
    cfg = SensorCfg(acc_bias=acc_bias, acc_noise=0.0071, gyro_bias=gyro_bias,
                    gyro_noise=0.0007, gyro_instab=2.0e-5,
                    tilt_deg=2.0, wobble_deg=wobble_deg, wobble_hz=1.2,
                    dt=dt, profile="cycloid")
    tr = ES.gen_truth_2(cfg, D, T_mot, seed=seed)
    down = max(0, tr["n_pre"] - int(press_early / dt))
    up = min(tr["N"], tr["n_pre"] + tr["n_mot"] + int(release_late / dt))
    return tr, down, up


def eval_case(label, D, T_mot, wobble, seeds=12, quiet=0.0, **kw):
    """quiet：把两端的静止段采样降噪？不需要 —— 仿真里两端本来就是纯静止。

    D=0（静止窗口）时相对误差没有定义，这时就报**绝对距离**：v1 会给出多少厘米。
    这正是用户报的那个毛病的量化证据，不能让它在打印层退化成 nan。
    """
    errs_v2, errs_v1, refusals, dists2, dists1 = [], [], [], [], []
    n_sens = 0
    for s in range(seeds):
        tr, dn, up = make(D, T_mot, wobble, seed=s, **kw)
        r = M2.measure(tr["acc"], tr["gyr"], dn, up, dt=DT)
        d1 = v1_distance(tr, dn, up, DT)  # v1 从不拒答，永远给个数
        dists1.append(d1)
        if not r["ok"]:
            refusals.append(r["reason"])
            continue
        dists2.append(r["distanceM"])
        if r.get("modelSensitive"):
            n_sens += 1
        if D > 0:
            errs_v2.append((r["distanceM"] - D) / D * 100)
            errs_v1.append((d1 - D) / D * 100)
    out = {"label": label, "D": D, "n_ok": len(dists2), "n_ref": len(refusals),
           "n_sens": n_sens,
           "median": statistics.median(dists2) * 100 if dists2 else float("nan"),
           "lo": min(dists2) * 100 if dists2 else float("nan"),
           "hi": max(dists2) * 100 if dists2 else float("nan"),
           "med2_abs": statistics.median(dists2) * 100 if dists2 else float("nan"),
           "med1_abs": statistics.median(dists1) * 100 if dists1 else float("nan"),
           "err2": statistics.median(errs_v2) if errs_v2 else float("nan"),
           "err2max": max(abs(e) for e in errs_v2) if errs_v2 else float("nan"),
           "err1": statistics.median(errs_v1) if errs_v1 else float("nan"),
           "err1max": max(abs(e) for e in errs_v1) if errs_v1 else float("nan"),
           "ref": refusals[0] if refusals else ""}
    return out


def row(r):
    mark = "★" if r["D"] == 0 else " "
    # 「模型敏感」那一列：仿真里注入了**完美常量**零偏，线性去趋势是精确的，
    # 所以这里每一次触发都是**误报**。把它显示出来，是为了别把这个判据吹得比实际好。
    sens = f" 模型敏感{r['n_sens']}/{r['n_ok']}"
    if r["D"] == 0:
        # 相对误差无定义 → 报绝对厘米数（真值是 0，所以"报了多少"就是"错了多少"）
        body = (f"v2 报 {r['med2_abs']:6.2f}cm          "
                f"| v1 报 {r['med1_abs']:7.2f}cm          "
                f"| v2 答 {r['n_ok']:2d}/拒 {r['n_ref']:2d}")
    else:
        body = (f"v2 中位 {r['err2']:+6.2f}% 最大 {r['err2max']:5.2f}%  "
                f"| v1 中位 {r['err1']:+7.2f}% 最大 {r['err1max']:6.2f}%  "
                f"| v2 答 {r['n_ok']:2d}/拒 {r['n_ref']:2d}")
    return f"{mark}{r['label']:>26}{r['D']*100:7.1f}cm  " + body + sens \
        + (f"  〔{r['ref']}〕" if r["ref"] else "")


def main():
    print("=" * 150)
    print("v2 vs v1（12 个种子；误差相对真值）")
    print("-" * 150)

    print("A · 静止窗口（用户报的毛病：按着不动也出距离）")
    for D, T, w in [(0.0, 4.0, 0.0), (0.0, 4.0, 2.0), (0.0, 8.0, 0.0)]:
        print(row(eval_case(f"D={D*100:.0f}cm T={T:.0f}s 转{w:.0f}°", D, T, w)))
    print()

    print("B · 正常滑动（按钮窗比运动长 0.55 s，两端各留静止）")
    for D, T in [(0.15, 1.5), (0.15, 0.8), (0.34, 1.5), (0.80, 1.5)]:
        print(row(eval_case(f"D={D*100:.0f}cm T={T:.1f}s", D, T, 1.0)))
    print()

    print("C · 运动中缓慢旋转（用户报的「对姿态要求也很严格」）")
    for w in (0.0, 2.0, 5.0, 10.0, 20.0):
        print(row(eval_case(f"D=15cm T=1.5s 转{w:.0f}°", 0.15, 1.5, w)))
    print()

    print("D · 陀螺零偏未校准（静止段极短 / 温漂）")
    for gb in (0.0, 0.001, 0.003):
        print(row(eval_case(f"gyro_bias={gb:.3f}", 0.15, 1.5, 2.0, gyro_bias=gb)))
    print()

    print("E · 加速度计零偏（去重力只能扣掉一部分）")
    for ab in (0.0, 0.092, 0.5):
        print(row(eval_case(f"acc_bias={ab:.3f}", 0.15, 1.5, 1.0, acc_bias=ab)))
    print()

    print("F · 只隔离姿态参考系那一项（真值窗口 + 无条件去趋势）")
    print(f"{'工况':>26}{'v1姿态(体参考)':>18}{'v2姿态(实测参考)':>18}")
    for w in (0.0, 2.0, 5.0, 10.0, 20.0):
        e1s, e2s = [], []
        for s in range(12):
            tr, dn, up = make(0.15, 1.5, w, seed=s)
            e1s.append((v1_attitude_off(tr, tr["n_pre"], tr["n_pre"] + tr["n_mot"], DT) / 0.15 - 1) * 100)
            r = M2.measure(tr["acc"], tr["gyr"], dn, up, dt=DT)
            if r["ok"]:
                e2s.append((r["distanceM"] - 0.15) / 0.15 * 100)
        print(f"{f'转{w:.0f}° 中位误差':>26}{statistics.median(e1s):17.2f}%{statistics.median(e2s):17.2f}%")


if __name__ == "__main__":
    main()
