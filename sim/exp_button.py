"""按钮手动边界（按住=开始，松手=结束）的可行性验证

用户的方案：不用算法去"检测"运动起止，改用一个大按钮，用户按住即开始、松手即结束，
DOWN/UP 事件的时刻直接当作积分区间 [lo, hi]。

要验证的核心命题（这句话成立与否决定方案成败）：

    积分窗口的边界落在「静止段」内部时，误差 ≈ 0。
    因为静止段的线性加速度 ≈ 0，多积一段不产生位移。

若成立 → 按钮不需要"精确"，只需要"落在静止段里"；
      → 用户天生的反应时间（看到按钮响应才推、推到位才松手）方向恰好是安全的；
      → 这比做亚采样变化点检测更鲁棒，且两端都带静止段会让 ZUPT 更干净。

要打掉的怀疑：
  (1) 窗口两端带静止段 → PCA 主轴会不会被静止段噪声拉偏？
  (2) 窗口变长 → ZUPT 的线性去趋势斜率变了会不会引入偏差？
  (3) 按钮时间戳抖动多大才会开始咬人？
  (4) 与"亚采样检测器"相比谁更好？

用 K40 实测参数（acc_noise 0.0071 / gyro_bias 0.0007°/s / dt 9.513ms）。
"""

import math
import random
import statistics
import sys

import pipeline as P
from pipeline import DT, G, D_TRUE, Params, SensorCfg, gen_truth

# ---- K40 实测参数（docs/health-a.md）------------------------------------
# profile 用 cycloid（摆线速度）：位移 ∝ (1-cos)，加速度首尾为 0 —— 真实手推的形态。
# 对照：sine 剖面的加速度首尾是 ±a_peak（不连续），会让梯形积分在跳变处多算半步，
#       那是"真值生成的不物理"，不是算法问题。见 exp_diag 的 6c。
K40 = dict(acc_bias=0.092, acc_noise=0.0071, gyro_bias=math.radians(0.0007),
           gyro_noise=0.0007, gyro_instab=2.0e-5, acc_scale=1.002, dt=0.009513,
           profile="cycloid")

# ---- 一台"标准设备"（教科书参数），用来对照"好机器 vs 差机器"--------------
BOOK = dict(acc_bias=0.05, acc_noise=0.02, gyro_bias=math.radians(0.5),
            gyro_noise=0.002, gyro_instab=3.5e-5, acc_scale=1.000, dt=0.01,
            profile="cycloid")


def rel(v):
    return (v - D_TRUE) / D_TRUE * 100.0


def run_manual(truth, p, lo, hi, pca_lo=None, pca_hi=None, qs=None, zupt=True):
    """按钮窗口版管线：积分区间由外部给定（lo, hi），其余与 run_pipeline 一致。

    pca_lo/pca_hi 缺省与积分区间相同；若分开给，可以单独考察"PCA 用运动段、积分用按钮窗"。
    qs 给定非 None 时用外部姿态（诊断用，比如真值姿态）；zupt=False 跳过线性去趋势。
    """
    dt = p.assumed_dt
    if qs is None:
        prof = P.build_profile(truth, p)
        qs, _ = P.solve_attitude(truth, prof.gyro_bias, prof.gravity_ref, p, dt=dt)
    a_lin = P.remove_gravity(truth, qs, g_mag=G)

    plo = lo if pca_lo is None else pca_lo
    phi = hi if pca_hi is None else pca_hi
    u, lam1, lam2 = P.pca_axis(a_lin[plo:phi + 1])
    axis = u
    if P.vdot(u, (1.0, 0.0, 0.0)) < 0:
        axis = P.vscale(axis, -1.0)
    a1 = [P.vdot(a_lin[i], axis) for i in range(len(a_lin))]

    v_raw, p_raw = P.integrate(a1, lo, hi, dt=dt)
    if not zupt:
        return abs(p_raw[-1]), (lam1 / lam2 if lam2 > 0 else float("inf")), abs(v_raw[-1])
    T = (hi - lo) * dt
    da = v_raw[-1] / T if T > 1e-9 else 0.0
    a1c = [a1[i] - da if lo <= i <= hi else a1[i] for i in range(len(a1))]
    _, p_z = P.integrate(a1c, lo, hi, dt=dt)
    return abs(p_z[-1]), (lam1 / lam2 if lam2 > 0 else float("inf")), abs(v_raw[-1])


def samples(sec, dt):
    return int(round(sec / dt))


# ==========================================================================
# 实验 1：起点锚点敏感度（终点固定在真值终点）
# ==========================================================================
def exp_start():
    print("=" * 96)
    print("1 · 起点锚点：把积分起点从真实起点挪动，看距离误差（终点固定在真值终点）")
    print("=" * 96)
    cfg = SensorCfg(**K40)
    p = Params(static_mode="adaptive", assumed_dt=cfg.dt)

    # 用 8 个种子取中位，压掉单次随机
    print(f"  {'起点相对真值':<18}{'提前/滞后':>12}{'距离中位':>11}{'误差中位':>11}{'最差':>10}")
    for k in [-samples(1.0, cfg.dt), -samples(0.5, cfg.dt), -samples(0.2, cfg.dt),
              -samples(0.05, cfg.dt), 0, samples(0.05, cfg.dt), samples(0.1, cfg.dt),
              samples(0.2, cfg.dt), samples(0.5, cfg.dt)]:
        errs = []
        for s in range(8):
            t = gen_truth(cfg, seed=3000 + s)
            lo = t["n_pre"] + k
            hi = t["n_pre"] + t["n_mot"]
            d, _, _ = run_manual(t, p, lo, hi)
            errs.append(rel(d))
        errs.sort()
        tag = f"{k * cfg.dt * 1000:+.0f} ms"
        side = "早（安全区）" if k < 0 else ("晚（切进运动段）" if k > 0 else "正好")
        print(f"  {tag:<18}{side:>12}{D_TRUE*(1+statistics.median(errs)/100)*100:>10.2f}cm"
              f"{statistics.median(errs):>10.2f}%{errs[-1]:>9.2f}%")
    print("\n  → 关键：'早于'曲线是平的（≈0），'晚于'曲线陡峭。")
    print()


# ==========================================================================
# 实验 2：终点锚点敏感度（起点固定在真值起点）
# ==========================================================================
def exp_end():
    print("=" * 96)
    print("2 · 终点锚点：把积分终点从真实终点挪动（起点固定在真值起点）")
    print("=" * 96)
    cfg = SensorCfg(**K40)
    p = Params(static_mode="adaptive", assumed_dt=cfg.dt)

    print(f"  {'终点相对真值':<18}{'提前/滞后':>12}{'距离中位':>11}{'误差中位':>11}{'最差':>10}")
    for k in [-samples(0.5, cfg.dt), -samples(0.2, cfg.dt), -samples(0.05, cfg.dt), 0,
              samples(0.05, cfg.dt), samples(0.1, cfg.dt), samples(0.2, cfg.dt),
              samples(0.5, cfg.dt), samples(1.0, cfg.dt)]:
        errs = []
        for s in range(8):
            t = gen_truth(cfg, seed=3000 + s)
            lo = t["n_pre"]
            hi = t["n_pre"] + t["n_mot"] + k
            d, _, _ = run_manual(t, p, lo, hi)
            errs.append(rel(d))
        errs.sort()
        tag = f"{k * cfg.dt * 1000:+.0f} ms"
        side = "晚（安全区）" if k > 0 else ("早（切进运动段）" if k < 0 else "正好")
        print(f"  {tag:<18}{side:>12}{D_TRUE*(1+statistics.median(errs)/100)*100:>10.2f}cm"
              f"{statistics.median(errs):>10.2f}%{errs[-1]:>9.2f}%")
    print("\n  → 端点的不对称性与起点镜像：'晚于'安全，'早于'危险。")
    print()


# ==========================================================================
# 实验 3：按钮时间戳抖动 —— 需要多大的提前量才安全
# ==========================================================================
def exp_jitter():
    print("=" * 96)
    print("3 · 按钮时间戳抖动：真实操作中按钮时刻会抖，看需要多少'安全余量'")
    print("=" * 96)
    cfg = SensorCfg(**K40)
    p = Params(static_mode="adaptive", assumed_dt=cfg.dt)

    print("  假设用户按下比真实起点早 margin、松手比真实终点晚同样的 margin，")
    print("  但两者都叠加一个 N(0, σ) 的抖动。抖到超出 margin 才切进运动段。")
    print()
    print(f"  {'安全余量':<12}" + "".join(f"{'σ=' + str(s) + 'ms':>13}" for s in [0, 10, 20, 30, 50]))
    for margin_ms in [0, 20, 50, 100, 200, 300]:
        row = ""
        for sig_ms in [0, 10, 20, 30, 50]:
            errs = []
            for s in range(12):
                rng = random.Random(4000 + s)
                t = gen_truth(cfg, seed=4000 + s, rng=rng)
                m = samples(margin_ms / 1000.0, cfg.dt)
                js = samples(sig_ms / 1000.0, cfg.dt)
                kj = int(round(rng.gauss(0, js))) if js > 0 else 0
                jk = int(round(rng.gauss(0, js))) if js > 0 else 0
                lo = t["n_pre"] - m + kj
                hi = t["n_pre"] + t["n_mot"] + m + jk
                d, _, _ = run_manual(t, p, lo, hi)
                errs.append(rel(d))
            errs.sort()
            row += f"{statistics.median(errs):>12.2f}%"
        print(f"  {str(margin_ms) + ' ms':<12}{row}")
    print("\n  → 只要安全余量 ≳ 3σ，抖动完全被吸收。人的反应时间天然提供 100ms 级余量。")
    print()


# ==========================================================================
# 实验 4：两端都带静止段，会不会污染 PCA / ZUPT？
# ==========================================================================
def exp_side_effects():
    print("=" * 96)
    print("4 · 副作用体检：窗口两端带静止段，会不会把 PCA 主轴或 ZUPT 去趋势带偏？")
    print("=" * 96)
    cfg = SensorCfg(**K40)
    p = Params(static_mode="adaptive", assumed_dt=cfg.dt)

    print(f"  {'余量':<10}{'积分窗=PCA窗':>16}{'误差中位':>11}{'λ1/λ2':>10}   "
          f"{'PCA用真值运动段':>18}{'误差中位':>11}{'λ1/λ2':>10}")
    for margin_ms in [0, 50, 100, 200, 400, 800]:
        a_errs, a_lam, b_errs, b_lam = [], [], [], []
        m = samples(margin_ms / 1000.0, cfg.dt)
        for s in range(10):
            t = gen_truth(cfg, seed=5000 + s)
            lo, hi = t["n_pre"] - m, t["n_pre"] + t["n_mot"] + m
            d, lam, _ = run_manual(t, p, lo, hi)
            a_errs.append(rel(d)); a_lam.append(lam)
            # PCA 改用真值运动段
            d2, lam2, _ = run_manual(t, p, lo, hi,
                                     pca_lo=t["n_pre"], pca_hi=t["n_pre"] + t["n_mot"])
            b_errs.append(rel(d2)); b_lam.append(lam2)
        a_errs.sort(); b_errs.sort()
        print(f"  {str(margin_ms) + ' ms':<10}{statistics.median(a_errs):>15.2f}%"
              f"{statistics.median(a_lam):>10.0f}   "
              f"{statistics.median(b_errs):>17.2f}%{statistics.median(b_lam):>10.0f}")
    print("\n  → 若两列接近，说明静止段对 PCA 无害（运动段能量远大于静止段噪声）。")
    print()


# ==========================================================================
# 实验 5：按钮方案 vs 算法检测方案（同一批数据）
# ==========================================================================
def exp_vs_detector():
    print("=" * 96)
    print("5 · 对决：按钮手动边界 vs 亚采样前的滑窗检测器（同一批合成数据）")
    print("=" * 96)
    print(f"  {'设备':<26}{'滑窗检测(现方案)':>18}{'按钮(安全余量200ms)':>21}{'按钮(余量0)':>14}")

    for name, kw, margin_ms in [("K40 实测（摆线剖面）", {**K40, "profile": "cycloid"}, 200),
                                ("教科书（摆线剖面）", {**BOOK, "profile": "cycloid"}, 200)]:
        cfg = SensorCfg(**kw)
        p_auto = Params(static_mode="adaptive", assumed_dt=cfg.dt)
        p_btn = Params(static_mode="adaptive", assumed_dt=cfg.dt)

        auto_errs, btn200_errs, btn0_errs = [], [], []
        for s in range(20):
            t = gen_truth(cfg, seed=6000 + s)
            r = P.run_pipeline(t, cfg, p_auto, seed=6000 + s)
            if r.get("ok", True):
                auto_errs.append(rel(r["dist"]))
            m = samples(margin_ms / 1000.0, cfg.dt)
            d, _, _ = run_manual(t, p_btn, t["n_pre"] - m, t["n_pre"] + t["n_mot"] + m)
            btn200_errs.append(rel(d))
            d0, _, _ = run_manual(t, p_btn, t["n_pre"], t["n_pre"] + t["n_mot"])
            btn0_errs.append(rel(d0))
        auto_errs.sort(); btn200_errs.sort(); btn0_errs.sort()
        a_med = statistics.median(auto_errs) if auto_errs else float("nan")
        a_max = auto_errs[-1] if auto_errs else float("nan")
        print(f"  {name:<26}{f'{a_med:.2f}% / {a_max:.2f}%':>18}"
              f"{f'{statistics.median(btn200_errs):.2f}% / {btn200_errs[-1]:.2f}%':>21}"
              f"{statistics.median(btn0_errs):>13.2f}%")
    print("\n  格式：中位 / 最大。")
    print()


# ==========================================================================
# 实验 6：诊断 —— "留余量后 +1~2% 的正偏差"到底是谁带来的
# ==========================================================================
def exp_diag():
    print("=" * 96)
    print("6a · 种子敏感性：同一物理配置换种子集，中位数会不会漂（判断实验1/3/4 的差异真假）")
    print("=" * 96)
    cfg = SensorCfg(**K40)
    p = Params(static_mode="adaptive", assumed_dt=cfg.dt)
    print(f"  {'种子基':<10}" + "".join(f"{'余量' + str(m) + 'ms':>14}" for m in [0, 100, 200, 400]))
    for seed0 in (3000, 4000, 5000, 7000, 9000):
        row = ""
        for margin_ms in (0, 100, 200, 400):
            errs = []
            for s in range(40):
                t = gen_truth(cfg, seed=seed0 + s)
                m = samples(margin_ms / 1000.0, cfg.dt)
                d, _, _ = run_manual(t, p, t["n_pre"] - m, t["n_pre"] + t["n_mot"] + m)
                errs.append(rel(d))
            errs.sort()
            row += f"{statistics.median(errs):>13.3f}%"
        print(f"  {str(seed0):<10}{row}")
    print()

    print("=" * 96)
    print("6b · 来源隔离：逐层换成「理想件」，看正偏差消失在谁身上")
    print("=" * 96)

    def dist_from_a1(a1, lo, hi, dt, zupt=True):
        v, pp = P.integrate(a1, lo, hi, dt=dt)
        if not zupt:
            return abs(pp[-1])
        T = (hi - lo) * dt
        da = v[-1] / T if T > 1e-9 else 0.0
        a1c = [a1[i] - da if lo <= i <= hi else a1[i] for i in range(len(a1))]
        _, pz = P.integrate(a1c, lo, hi, dt=dt)
        return abs(pz[-1])

    print(f"  {'余量':<9}{'①全管线':>11}{'②真值姿态':>11}{'③不去趋势':>11}"
          f"{'④真姿态+不去势':>16}{'⑤纯真值积分':>14}")
    for margin_ms in (0, 100, 200, 400, 800):
        arms = [[], [], [], [], []]
        m = samples(margin_ms / 1000.0, cfg.dt)
        for s in range(30):
            t = gen_truth(cfg, seed=7000 + s)
            lo, hi = t["n_pre"] - m, t["n_pre"] + t["n_mot"] + m
            arms[0].append(rel(run_manual(t, p, lo, hi)[0]))
            arms[1].append(rel(run_manual(t, p, lo, hi, qs=t["q_true"])[0]))
            arms[2].append(rel(run_manual(t, p, lo, hi, zupt=False)[0]))
            arms[3].append(rel(run_manual(t, p, lo, hi, qs=t["q_true"], zupt=False)[0]))
            arms[4].append(rel(dist_from_a1(t["a_true"], lo, hi, cfg.dt)))
        med = [statistics.median(a) for a in arms]
        print(f"  {str(margin_ms) + ' ms':<9}" + "".join(f"{v:>10.3f}%" for v in med[:3])
              + f"{med[3]:>15.3f}%{med[4]:>13.3f}%")
    print()
    print("  ① 估计姿态 + ZUPT（现状）    ② 换真值姿态（隔离姿态误差）")
    print("  ③ 估计姿态但跳过线性去趋势    ④ 真值姿态且不去趋势")
    print("  ⑤ 直接用真值加速度序列积分（不含任何传感器误差，只含离散化）")
    print()
    print("  判读：② ≈ ① → 姿态误差不是原因；③ 更小 → 去趋势贡献一点；")
    print("        ⑤ 若也随余量上升 → 偏差来自真值序列本身，与按钮方案无关（见 6c）。")
    print()

    print("=" * 96)
    print("6c · 决定性对照：运动剖面「加速度首尾是否连续」对窗口余量的敏感性")
    print("=" * 96)
    print("  sine   : 正弦速度剖面，a(0)=+a_peak —— 加速度在起点不连续（不物理）")
    print("  cycloid: 摆线速度剖面，a(0)=a(T)=0 —— 真实推手的形态")
    print()
    print(f"  {'余量':<10}{'sine：误差中位':>17}{'cycloid：误差中位':>20}{'cycloid 最大':>14}")
    for margin_ms in (0, 100, 200, 400, 800):
        cells = []
        for prof in ("sine", "cycloid"):
            errs = []
            for s in range(30):
                t = gen_truth(SensorCfg(**{**K40, "profile": prof}), seed=7000 + s)
                m = samples(margin_ms / 1000.0, cfg.dt)
                d, _, _ = run_manual(t, p, t["n_pre"] - m, t["n_pre"] + t["n_mot"] + m)
                errs.append(rel(d))
            errs.sort()
            cells.append((statistics.median(errs), errs[-1]))
        print(f"  {str(margin_ms) + ' ms':<10}{cells[0][0]:>16.3f}%{cells[1][0]:>19.3f}%"
              f"{cells[1][1]:>13.3f}%")
    print()
    print("  → cycloid 下误差与余量无关（平的）；sine 下的上升是真值不连续造的伪影。")
    print("     真实手推的加速度首尾必然连续，所以按钮方案的窗口余量是真正安全的。")
    print()


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "1"):
        exp_start()
    if which in ("all", "2"):
        exp_end()
    if which in ("all", "3"):
        exp_jitter()
    if which in ("all", "4"):
        exp_side_effects()
    if which in ("all", "5"):
        exp_vs_detector()
    if which in ("all", "6", "diag"):
        exp_diag()


if __name__ == "__main__":
    main()
