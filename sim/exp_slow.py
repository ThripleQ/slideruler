"""慢速短距滑动：算法到底还有没有分辨力？

背景：20261011_090006 真机采集（K40），用户报真值 15 cm，按住 4.4 s。
算法给 12.24 cm（−18%）。用 K40 参数在仿真里复现这个工况，看差距是
"算法在慢速区就是不行" 还是 "这次采集有别的毛病"。

用法：python exp_slow.py
"""

import math
import random
import statistics

import pipeline as P
from pipeline import (DT, G, DeviceProfile, Params, SensorCfg, q_from_euler,
                      qconj, qmul, qrotT, vdot, vnorm, vscale, vunit)
from exp_button import K40, samples


def gen_truth_2(cfg, D, T_mot, seed=0):
    """和 pipeline.gen_truth 同构，但把 D_TRUE / MOTION_S 参数化。"""
    rng = random.Random(seed)
    dt = cfg.dt
    n_pre = int(round(2.0 / dt))
    n_mot = int(round(T_mot / dt))
    n_post = int(round(1.2 / dt))
    N = n_pre + n_mot + n_post

    def _s(tau):
        return D * (tau / T_mot - math.sin(2 * math.pi * tau / T_mot) / (2 * math.pi))

    def _a(tau):
        return (2.0 * D / T_mot) * (math.pi / T_mot) * math.sin(2 * math.pi * tau / T_mot)

    s = [0.0] * (N + 1)
    a_true = [0.0] * (N + 1)
    for i in range(N + 1):
        tau = (i - n_pre) * dt
        if n_pre <= i <= n_pre + n_mot:
            a_true[i] = _a(tau)
        s[i] = 0.0 if i <= n_pre else (D if i >= n_pre + n_mot else _s(tau))

    tilt = math.radians(cfg.tilt_deg)
    amp = math.radians(cfg.wobble_deg)
    ww = 2 * math.pi * cfg.wobble_hz
    q_true = []
    for i in range(N + 1):
        t = i * dt
        u = (i - n_pre) * dt / T_mot
        env = math.sin(math.pi * u) ** 2 if 0.0 < u < 1.0 else 0.0
        k = amp * env
        roll = tilt + k * math.sin(ww * t)
        pitch = tilt * 0.6 + k * 0.7 * math.cos(ww * t * 0.8)
        q_true.append(q_from_euler(roll, pitch, 0.0))

    gyro_clean = []
    for i in range(N + 1):
        j = min(i + 1, N)
        dq = qmul(qconj(q_true[i]), q_true[j])
        gyro_clean.append((2.0 * dq[1] / dt, 2.0 * dq[2] / dt, 2.0 * dq[3] / dt))

    acc, gyr = [], []
    instability = [0.0, 0.0, 0.0]
    for i in range(N + 1):
        for a in range(3):
            instability[a] += cfg.gyro_instab * math.sqrt(dt) * rng.gauss(0, 1)
        f_world = (a_true[i], 0.0, G)
        f_body = qrotT(q_true[i], f_world)
        acc.append(tuple(cfg.acc_scale * f_body[a] + cfg.acc_bias
                         + cfg.acc_noise * rng.gauss(0, 1) for a in range(3)))
        gyr.append(tuple(gyro_clean[i][a] + cfg.gyro_bias + instability[a]
                         + cfg.gyro_noise * rng.gauss(0, 1) for a in range(3)))
    return {"N": N, "n_pre": n_pre, "n_mot": n_mot, "s": s, "a_true": a_true,
            "q_true": q_true, "acc": acc, "gyr": gyr, "dt": dt}


def make_cfg(dt):
    return SensorCfg(acc_bias=K40["acc_bias"], acc_noise=K40["acc_noise"],
                     gyro_bias=K40["gyro_bias"], gyro_noise=0.0007,
                     gyro_instab=2.0e-5, tilt_deg=2.0, wobble_deg=1.0,
                     wobble_hz=1.2, dt=dt, profile="cycloid")


def pipeline_distance(tr, dt, margin_s=0.3, g_mag=None, cond_zupt=True):
    """当前算法的完整一维距离（与 slide_pipeline.py 同构）。

    g_mag=None → 用本机静止段实测的重力模长（正确做法）。
    g_mag=G    → 复现"写死 9.81"的旧行为（对照用）。
    cond_zupt=False → 复现"无条件 ZUPT"的旧行为（对照用）。
    """
    n = tr["N"]
    prof = P.build_profile(tr)
    p = Params(static_mode="adaptive", assumed_dt=dt)
    qs, _ = P.solve_attitude(tr, prof.gyro_bias, prof.gravity_ref, p, dt=dt)
    a_lin = P.remove_gravity(tr, qs, g_mag=prof.gravity_mag if g_mag is None else g_mag)

    lo = max(0, tr["n_pre"] - int(margin_s / dt))
    hi = min(n, tr["n_pre"] + tr["n_mot"] + int(margin_s / dt))

    u, lam1, lam2 = P.pca_axis(a_lin[lo:hi + 1])
    axis = u if vdot(u, (1.0, 0.0, 0.0)) >= 0 else vscale(u, -1.0)
    a1 = [vdot(a_lin[i], axis) for i in range(n + 1)]

    cal = range(max(0, lo - int(2.0 / dt)), lo)
    sigma = statistics.pstdev([a1[i] for i in cal])
    b = abs(statistics.fmean(a1[i] for i in cal))

    if cond_zupt:
        d, _, _, _ = P.integrate_cond_zupt(a1, lo, hi, dt, bias_bound=b, sigma=sigma)
    else:
        _, pos = P.integrate(a1, lo, hi, dt=dt)
        T = (hi - lo) * dt
        v, _ = P.integrate(a1, lo, hi, dt=dt)
        da = v[-1] / T if T > 1e-9 else 0.0
        a1c = [a1[i] - da if lo <= i <= hi else a1[i] for i in range(n + 1)]
        _, pz = P.integrate(a1c, lo, hi, dt=dt)
        d = abs(pz[-1])
    return d


def one(D, T, dt, margin_s=0.3, seed=0):
    cfg = make_cfg(dt)
    tr = gen_truth_2(cfg, D, T, seed=seed)
    d = pipeline_distance(tr, dt, margin_s=margin_s)
    return d / D - 1.0, d


def main():
    dt = 0.0095136
    print("K40 参数 · 按钮窗 +0.3 s 余量 · 摆线剖面（a 首尾为 0）")
    print("每格 = 20 次随机种子的误差中位数 / 最大 |误差|")
    print()
    print(f"{'位移\\时长':>10}" + "".join(f"{f'{t:.1f}s':>16}" for t in (0.8, 1.5, 2.3, 3.5, 4.5)))
    for D in (0.15, 0.34, 0.80):
        row = []
        for T in (0.8, 1.5, 2.3, 3.5, 4.5):
            errs = []
            for seed in range(20):
                e, _ = one(D, T, dt, seed=seed)
                errs.append(e)
            med = statistics.median(errs)
            mx = max(abs(e) for e in errs)
            row.append(f"{med*100:+6.1f}/{mx*100:5.1f}")
        print(f"{D*100:7.0f}cm " + "".join(f"{c:>16}" for c in row))
    print()
    print("（每格：中位误差 % / 最大绝对误差 %。摆线剖面的 a_peak = 2πD/T²·... 已含在内）")
    print()

    # 15 cm / 4.4 s：两个 bug 各自的贡献（各 20 种子，中位 / 最大 |误差|）
    print("本次工况复现：15 cm / 4.4 s —— 两个缺陷各自的贡献")
    print(f"{'配置':>34}{'误差中位':>12}{'最大':>10}")

    def sweep(label, g_mag, cond):
        es = []
        for seed in range(20):
            cfg = make_cfg(dt)
            tr = gen_truth_2(cfg, 0.15, 4.4, seed=seed)
            d = pipeline_distance(tr, dt, g_mag=g_mag, cond_zupt=cond)
            es.append(abs(d / 0.15 - 1) * 100)
        es.sort()
        print(f"{label:>34}{statistics.median(es):11.2f}%{es[-1]:9.2f}%")

    sweep("旧：g_mag=9.81 + 无条件 ZUPT", G, False)
    sweep("只修 g_mag（仍无条件 ZUPT）", None, False)
    sweep("只改条件 ZUPT（仍 g_mag=9.81）", G, True)
    sweep("当前：g_mag=‖g_meas‖ + 条件 ZUPT", None, True)

    # 不同时长的残余误差（当前算法）
    print()
    print("当前算法：随按住时长增长的残余误差（15 cm）")
    print(f"{'时长':>8}{'误差中位':>12}{'最大':>10}{'裸积分/去趋势各占多少':>24}")
    for T in (0.8, 1.2, 1.5, 2.3, 3.5, 4.4):
        es = sorted(abs(one(0.15, T, dt, seed=s)[0]) * 100 for s in range(20))
        print(f"{T:7.1f}s{statistics.median(es):11.2f}%{es[-1]:9.2f}%")


if __name__ == "__main__":
    main()
