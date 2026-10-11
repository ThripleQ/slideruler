"""slideruler 合成数据仿真总入口。

python run.py        跑全部实验
python run.py quick  只跑关键结论
"""

import math
import statistics
import sys

from pipeline import (DT, D_TRUE, MOTION_S, Params, SensorCfg, gen_truth,
                      run_pipeline)

IDEAL = dict(acc_bias=0, acc_noise=0, gyro_bias=0, gyro_noise=0, gyro_instab=0,
             wobble_deg=0, tilt_deg=0)


def rel(d):
    return (d - D_TRUE) / D_TRUE * 100.0


def one(cfg, p, seed):
    return run_pipeline(gen_truth(cfg, seed=seed), cfg, p, seed)


# ---------------------------------------------------------------- 1. 标称诊断
def sec_nominal():
    print("=" * 104)
    print("1 · 标称工况端到端（0.8 m / 1.5 s / 100 Hz，传感器参数取指南 §1.3）")
    print("=" * 104)
    cfg = SensorCfg()
    truth = gen_truth(cfg, seed=7)
    a_pk = max(abs(x) for x in truth["a_true"])
    print(f"  正弦速度剖面：v_peak {math.pi * D_TRUE / (2 * MOTION_S):.3f} m/s，"
          f"a_peak {a_pk:.3f} m/s²  ← 指南 §5 写的是 ~1.1 m/s²")
    for mode, desc in [("oracle", "真值窗口（隔离检测器）"),
                       ("linear_accel", "‖a_lin‖ 判据（推荐）"),
                       ("guide_magnitude", "指南 §7 判据一")]:
        p = Params(static_mode=mode)
        r = run_pipeline(truth, cfg, p, seed=7)
        if not r.get("ok", True):
            print(f"  {desc:<22} 检测失败：{r['reason']}")
            continue
        print(f"  {desc:<22} 窗口[{r['lo']:>3},{r['hi']:>3}]  姿态误差 {r['att_err_motion']:.3f}°  "
              f"λ₁/λ₂ {r['lam_ratio']:>6.1f}  v_err {r['v_err']:+.4f} m/s  →  "
              f"{r['dist'] * 100:.2f} cm ({rel(r['dist']):+.2f}%)")
    print()


# ----------------------------------------------------- 2. 边界定位：真正的瓶颈
def sec_boundary():
    print("=" * 104)
    print("2 · 误差的真正来源：积分起点定位（理想传感器 + 真值姿态 + 理想主轴）")
    print("=" * 104)
    cfg = SensorCfg(**IDEAL)
    truth = gen_truth(cfg, seed=7)
    a1 = [truth["a_true"][i] for i in range(len(truth["a_true"]))]

    def integ(lo, hi, drift=True):
        v = p = 0.0
        vs = [0.0]
        for i in range(lo + 1, hi + 1):
            v += (a1[i - 1] + a1[i]) / 2 * DT
            vs.append(v)
        da = vs[-1] / ((hi - lo) * DT) if drift else 0.0
        v = p = 0.0
        for i in range(lo + 1, hi + 1):
            vn = v + (a1[i - 1] + a1[i]) / 2 * DT - da * DT
            p += (v + vn) / 2 * DT
            v = vn
        return p

    print("  （两种偏移都按管线做法做线性漂移扣除）")
    print(f"  {'起点晚':<16}{'距离':>10}{'误差':>10}      {'终点晚':<16}{'距离':>10}{'误差':>10}")
    for k in [0, 1, 2, 5, 10, 30]:
        dl = integ(200 + k, 350)
        dr = integ(200, 350 + k)
        print(f"  {str(k * 10) + ' ms':<16}{dl * 100:>10.2f}{rel(dl):>9.2f}%      "
              f"{str(k * 10) + ' ms':<16}{dr * 100:>10.2f}{rel(dr):>9.2f}%")
    print("\n  起点晚 10 ms 就 -1.7%，晚 50 ms 就 -8.3%，晚 300 ms 直接报废；")
    print("  终点晚 300 ms 只差 0.7%。")
    print("  原因：起点处速度刚从 0 抬起，把 v=0 钉晚一点，等于把真实速度整段丢掉；")
    print("        终点处速度已经回到 0 附近，多积一段静止段不产生位移。")
    print("  → 只有起点边界需要采样级精度，终点可以用粗判据。")
    print()


# ------------------------------------------------------------ 3. 检测器对比
def sec_detector():
    print("=" * 104)
    print("3 · 静止检测判据：起步形态决定它是可用还是废掉")
    print("=" * 104)
    for profile, desc in [("sine", "起步瞬间加速度即峰值（对检测器最友好）"),
                          ("cycloid", "起步加速度为 0、缓慢起推（更像真实手）")]:
        print(f"\n  ### {desc}")
        cfg = SensorCfg(profile=profile)
        truth = gen_truth(cfg, seed=7)
        r0 = run_pipeline(truth, cfg, Params(static_mode="oracle"), seed=7)
        a_pk = max(abs(x) for x in truth["a_true"])
        print(f"      加速度峰值 {a_pk:.3f} m/s²，起步 50 ms 后为 {truth['a_true'][205]:.3f} m/s²")
        for mode, mdesc in [("guide_magnitude", "指南 §7 判据一 |‖f‖-9.81|<0.2"),
                            ("linear_accel", "去重力后 ‖a_lin‖<0.12"),
                            ("accel_variance", "窗内加速度 σ<0.05"),
                            ("oracle", "真值窗口")]:
            p = Params(static_mode=mode)
            r = run_pipeline(truth, cfg, p, seed=7)
            if not r.get("ok", True):
                print(f"      {mdesc:<30} 检测失败（整段被判为静止）")
                continue
            lag = r["lo"] - 200
            print(f"      {mdesc:<30} 起点滞后 {lag * 10:>4} ms   距离 {r['dist'] * 100:>6.2f} cm"
                  f"   误差 {rel(r['dist']):>+7.2f}%")
    print()


# -------------------------------------------------------------- 4. 机制拆解
def sec_mechanism():
    print("=" * 104)
    print("4 · 机制拆解（窗口用真值，隔离检测器；各 25 次随机种子）")
    print("=" * 104)
    for label, mk in [("全管线", lambda c, q: None),
                      ("关掉 ZUPT", lambda c, q: setattr(q, "use_zupt", False)),
                      ("关掉一维投影 PCA", lambda c, q: setattr(q, "use_pca", False)),
                      ("关掉姿态解算", lambda c, q: setattr(q, "use_attitude", False)),
                      ("关掉三者", lambda c, q: (setattr(q, "use_attitude", False),
                                                setattr(q, "use_pca", False),
                                                setattr(q, "use_zupt", False)))]:
        errs = []
        for s in range(25):
            cfg, p = SensorCfg(), Params(static_mode="oracle")
            mk(cfg, p)
            r = run_pipeline(gen_truth(cfg, seed=1000 + s), cfg, p, seed=1000 + s)
            if r.get("ok", True):
                key = "dist_zupt" if (p.use_zupt and r.get("dist_zupt")) else "dist_raw"
                errs.append(abs(rel(r[key])))
        if errs:
            print(f"  {label:<20} 绝对误差中位 {statistics.median(errs):6.2f}%   "
                  f"最大 {max(errs):6.2f}%")
    print()


# ---------------------------------------------------------------- 5. 敏感性
def sec_sweep():
    print("=" * 104)
    print("5 · 敏感性扫描（窗口用真值；把指南给的数字逐个推坏）")
    print("=" * 104)

    def scan(title, setter, values, n=25):
        print(f"\n  ### {title}")
        for v in values:
            errs = []
            for s in range(n):
                cfg, p = SensorCfg(), Params(static_mode="oracle")
                setter(cfg, p, v)
                r = run_pipeline(gen_truth(cfg, seed=1000 + s), cfg, p, seed=1000 + s)
                if r.get("ok", True):
                    key = "dist_zupt" if (p.use_zupt and r.get("dist_zupt")) else "dist_raw"
                    errs.append(abs(rel(r[key])))
            errs.sort()
            print(f"      {str(v):<12} 中位 {statistics.median(errs):6.2f}%   "
                  f"P90 {errs[min(len(errs) - 1, int(0.9 * len(errs)))]:6.2f}%   "
                  f"最大 {errs[-1]:6.2f}%")

    scan("加速度计白噪声 σ m/s²（指南 §2 声称 0.02）",
         lambda c, q, v: setattr(c, "acc_noise", v), [0.01, 0.02, 0.04, 0.08, 0.16])
    scan("加速度计零偏 m/s²（指南 §1.3 用 0.05）",
         lambda c, q, v: setattr(c, "acc_bias", v), [0.0, 0.05, 0.1, 0.2, 0.5])
    scan("陀螺零偏 °/s（指南 §2 声称 0.5，校准段可消）",
         lambda c, q, v: setattr(c, "gyro_bias", math.radians(v)), [0.5, 1.0, 2.0, 4.0])
    scan("陀螺零偏不稳定性 rad/s/√s（指南未给此数，校准消不掉）",
         lambda c, q, v: setattr(c, "gyro_instab", v),
         [0.0, 3.5e-5, 7e-5, 1.4e-4, 3.5e-4])
    scan("互补滤波 α（指南 §4 建议 0.02~0.05）",
         lambda c, q, v: setattr(q, "alpha", v), [0.005, 0.02, 0.05, 0.1, 0.3])
    scan("手部转动幅度 °（指南未列入误差预算）",
         lambda c, q, v: setattr(c, "wobble_deg", v), [0.0, 1.0, 2.0, 4.0, 8.0])
    scan("静态安装倾角 °（指南 §3 说极斜也没关系）",
         lambda c, q, v: setattr(c, "tilt_deg", v), [0.0, 2.0, 10.0, 25.0, 40.0])


# ------------------------------------------------ 6. 跨设备：自校准值多少钱
def sec_dt_basis():
    print("=" * 104)
    print("6 · 时间基：采样间隔每台机都不一样（K40 实测 9.513 ms），写死 100 Hz 会错多少")
    print("=" * 104)
    print("  （近理想传感器 + 真值窗口，隔离掉其它误差，只看时间基）")
    print(f"  {'真实 dt':>10}{'真实速率':>11}{'写死 10 ms':>13}{'误差':>10}"
          f"{'用实测 dt':>13}{'误差':>10}")
    for dt in [0.005, 0.009513, 0.0125, 0.02]:
        cfg = SensorCfg(dt=dt, **IDEAL)
        truth = gen_truth(cfg, seed=7)
        r_bad = run_pipeline(truth, cfg, Params(static_mode="oracle", assumed_dt=0.01), seed=7)
        r_ok = run_pipeline(truth, cfg, Params(static_mode="oracle", assumed_dt=dt), seed=7)
        print(f"  {dt * 1000:>8.3f}ms{1.0 / dt:>9.1f}Hz"
              f"{r_bad['dist'] * 100:>11.2f}cm{rel(r_bad['dist']):>9.2f}%"
              f"{r_ok['dist'] * 100:>11.2f}cm{rel(r_ok['dist']):>9.2f}%")
    print("\n  位移 ∝ dt²（两次积分都乘 dt）→ dt 错 5% 距离就错 10%。")
    print("  这不是小数：K40 实测 105.1 Hz，若按请求的 100 Hz 积分，先天 +10%。")
    print("  → 必须用 event.timestamp 的实测间隔。本管线里就是 Params.assumed_dt。")
    print()


def sec_devices():
    print("=" * 104)
    print("7 · 跨设备：把参数写死 / 每会话自校准 / 再加单姿态标度修正（各 25 个种子）")
    print("=" * 104)
    devices = [
        ("教科书（指南 §1.3 假设）",
         dict(acc_bias=0.05, acc_noise=0.02, gyro_bias=math.radians(0.5),
              gyro_noise=0.002, gyro_instab=3.5e-5, acc_scale=1.000, dt=0.01)),
        ("K40 实测（本机验机结果）",
         dict(acc_bias=0.092, acc_noise=0.0071, gyro_bias=math.radians(0.0007),
              gyro_noise=0.0007, gyro_instab=2.0e-5, acc_scale=1.002, dt=0.01)),
        ("廉价机（噪声大 / 零偏大 / 标度 −2%）",
         dict(acc_bias=0.35, acc_noise=0.06, gyro_bias=math.radians(3.0),
              gyro_noise=0.006, gyro_instab=1.4e-4, acc_scale=0.98, dt=0.01)),
        ("老机（陀螺漂移大 / 标度 +2%）",
         dict(acc_bias=0.10, acc_noise=0.03, gyro_bias=math.radians(1.0),
              gyro_noise=0.003, gyro_instab=3.5e-4, acc_scale=1.02, dt=0.01)),
    ]
    cfgs = [
        ("A 写死", Params(static_mode="linear_accel", static_lin_thr=0.12,
                        use_scale_corr=False, assumed_dt=0.01)),
        ("B 自校准", Params(static_mode="adaptive", use_scale_corr=False, assumed_dt=0.01)),
        ("C B+单姿态标度", Params(static_mode="adaptive", use_scale_corr=True, assumed_dt=0.01)),
    ]

    print(f"  {'设备':<32}{'真实 k':>8}" + "".join(f"{n + ' 中位':>13}{'最大':>9}" for n, _ in cfgs))
    stats = {}
    for name, kw in devices:
        stats[name] = {}
        row = ""
        for key, p in cfgs:
            errs = []
            for s in range(25):
                cfg = SensorCfg(**kw)
                r = run_pipeline(gen_truth(cfg, seed=1000 + s), cfg, p, seed=1000 + s)
                if r.get("ok", True):
                    errs.append(abs(rel(r["dist"])))
            if not errs:
                stats[name][key] = None
                row += f"{'—':>13}{'—':>9}"
                continue
            errs.sort()
            stats[name][key] = (statistics.median(errs), errs[-1])
            row += f"{statistics.median(errs):>12.2f}%{errs[-1]:>8.2f}%"
        print(f"  {name:<32}{kw['acc_scale']:>8.3f}{row}")

    print("\n  自校准估出的 |g| 与 k̂ —— 注意它并不等于真实标度：")
    print(f"      {'设备':<32}{'真实 k':>8}{'‖g_meas‖':>11}{'k̂=‖g‖/9.81':>13}{'k̂ 相对真值的错':>16}")
    kerr = {}
    for name, kw in devices:
        cfg = SensorCfg(**kw)
        r = run_pipeline(gen_truth(cfg, seed=1000), cfg, Params(static_mode="adaptive"), seed=1000)
        prof = r["profile"]
        khat = prof.scale_hat()
        kerr[name] = (khat / kw["acc_scale"] - 1) * 100
        print(f"      {name:<32}{kw['acc_scale']:>8.3f}{prof.gravity_mag:>11.4f}"
              f"{khat:>13.5f}{kerr[name]:>15.2f}%")

    worst = max(devices, key=lambda d: stats[d[0]]["A 写死"][1])
    wn = worst[0]
    print("\n  结论 A→B（每会话自校准）：每台机都变好，而且把灾难性长尾砍掉了 ——")
    print(f"        最大误差降幅最大的是『{wn}』："
          f"{stats[wn]['A 写死'][1]:.2f}% → {stats[wn]['B 自校准'][1]:.2f}%；"
          f"中位 {stats[wn]['A 写死'][0]:.2f}% → {stats[wn]['B 自校准'][0]:.2f}%。")
    print("        收益几乎全部来自『静止检测阈值按本机噪声地板导出』，不是来自传感器零偏/标度。")

    print(f"\n  结论 C（再加单姿态标度修正）：**这是陷阱，不要用。**")
    print(f"        静止段读到的 ‖g‖ 混着『标度 k』和『零偏沿重力方向的投影 b·ĝ』，单姿态分不开。")
    print(f"        k̂ 相对真值的误差（右表最后一列，"
          f"{min(kerr.values()):+.2f}%~{max(kerr.values()):+.2f}%）普遍大于它想修的标度误差，")
    print("        于是 C 列成绩基本是随机的：有时碰巧抵消了检测器的正向偏差而显得好看，")
    print("        有时（廉价机）直接修坏。它不是校准，是一次带偏的扰动。")

    print("\n  两条附注：")
    print("        · 分离 k 与 b 的唯一办法是多姿态（≥4 个不共面朝向）球面拟合，见 fit_bias.py --selftest。")
    print("        · A/B 两列还有一个共同的正向偏差 ≈ +0.6%：滑窗检测器在滑动结束后仍多判约一个窗长，")
    print("          把终点边界推后了。这属于『边界检测』问题，不属于校准，留给下一步的亚采样检测器。")
    print()


def main():
    quick = len(sys.argv) > 1 and sys.argv[1] == "quick"
    sec_nominal()
    sec_boundary()
    sec_detector()
    sec_dt_basis()
    if not quick:
        sec_devices()
        sec_mechanism()
        sec_sweep()


if __name__ == "__main__":
    sys.exit(main())
