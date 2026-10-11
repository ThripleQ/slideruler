"""在物理一致性仿真上跑算法 —— 建立"可控真值"下的对照实验。

和真机的区别：这里的真值精确到 1e-9 m，且可以单独打开一种误差。
于是"这个改进到底有没有用"变成可回答的问题，而不是靠 3 条肉眼估的 15 cm 去猜。

每次只改变一个主要假设（开发纪律），所以开关全部显式：
    --v3       关掉 v4 的两项改动（杠杆臂 + 运动段禁用加速度修正），回到 v3
    --lever 0  只关杠杆臂
    --nogate   只关运动段禁用修正
    --mode button  运动区间改由按钮界定

用法：
    python sim/exp_sim.py --clean            # 无噪声，只验实现正确性
    python sim/exp_sim.py --seeds 30         # 有噪声，多种子统计
"""
import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from simgen import SimCfg, scenario, gen, down_i_of
import measure2


def run_case(spec, seed, mode, lever, gate, r_vec):
    spec["cfg"].seed = seed
    r = gen(spec["cfg"], spec["head_pos"], spec["att"], spec["r_body"], spec["t_end"])
    dt = 1.0 / spec["cfg"].hz
    dn = down_i_of(r["acc_ts"], spec["tA"] - 0.25)
    up = down_i_of(r["acc_ts"], spec["tB"] + 0.25)
    m = measure2.measure(r["acc"], r["gyr"], dn, up, dt=dt, merge_mode=mode,
                         head_r=r_vec if lever else None,
                         gate_motion=(gate if type(gate) is bool else True))
    return m


# (场景名, 参数, 真值 cm, 说明)
SCENARIOS = [
    ("static", dict(T=5.0), 0.0, "全程不动 —— 专测假距离"),
    ("jitter_only", dict(T=5.0, amp_deg=(2, 2, 2)), 0.0, "全程手抖，测量头不动"),
    ("static_rotate", dict(T=5.0, rot_deg=30), 0.0, "测量头不动、手机绕它转 30°"),
    ("slide", dict(T=5.0, d=0.15), 15.0, "理想直线滑动"),
    ("slide_rotate", dict(T=5.0, d=0.15, rot_deg=15), 15.0, "滑动 + 手机转 15°"),
    ("slide_wobble", dict(T=5.0, d=0.15, amp_deg=(3, 3, 4), drift_deg_per_s=(0, 0, 2),
                          only_during_slide=True), 15.0, "滑动时手抖（真机难度）"),
    ("slide_wobble", dict(T=5.0, d=0.15, amp_deg=(6, 6, 8), drift_deg_per_s=(0, 0, 4),
                          only_during_slide=True), 15.0, "滑动时手抖 ×2"),
    ("slide_wobble", dict(T=5.0, d=0.15, amp_deg=(3, 3, 4), drift_deg_per_s=(0, 0, 2)),
     15.0, "全程手抖（最坏情况）"),
]

CLEAN = dict(acc_noise=0.0, gyr_noise=0.0, acc_bias=(0, 0, 0), gyr_bias=(0, 0, 0),
             acc_scale=(0, 0, 0), acc_misalign=0.0)
# 噪声档：取 K40 / lsm6dso 实测量级
NOISY = dict(acc_noise=0.010, gyr_noise=0.002,
             acc_bias=(0.020, -0.010, 0.030), gyr_bias=(0.002, -0.001, 0.003),
             acc_scale=(0.011, 0.011, 0.011), acc_misalign=0.004)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clean", action="store_true")
    ap.add_argument("--seeds", type=int, default=1)
    ap.add_argument("--lever", type=float, default=0.09)
    ap.add_argument("--mode", default="press", choices=["press", "button", "overlap"])
    ap.add_argument("--v3", action="store_true", help="回到 v3：无杠杆臂 + 运动段仍做修正")
    ap.add_argument("--nolever", action="store_true")
    ap.add_argument("--nogate", action="store_true")
    a = ap.parse_args()

    lever = not (a.v3 or a.nolever or a.lever == 0)
    gate = not (a.v3 or a.nogate)
    base = CLEAN if a.clean else NOISY
    r_vec = (0.0, -a.lever, 0.0)

    print("杠杆臂 %s  运动段禁修正 %s  模式 %s  噪声 %s  种子 %d"
          % ("开(%.0fcm)" % (a.lever * 100) if lever else "关",
             "开" if gate else "关", a.mode, "无" if a.clean else "有", a.seeds))
    print("%-15s %6s %9s %8s %7s %6s  %s"
          % ("场景", "真值cm", "估cm", "误差cm", "相对", "拒答%", "说明"))
    print("-" * 100)

    for name, kw, truth_cm, desc in SCENARIOS:
        ests, rej = [], 0
        for s in range(a.seeds):
            spec = scenario(name, cfg=SimCfg(**base), r_body=r_vec, **kw)
            m = run_case(spec, s + 1, a.mode, lever, gate, r_vec)
            if m["ok"]:
                ests.append(m["distanceM"])
            else:
                rej += 1
        if ests:
            mu = statistics.fmean(ests)
            sd = statistics.pstdev(ests) if len(ests) > 1 else 0.0
            err = mu - truth_cm / 100.0
            print("%-15s %6.2f %9.3f %8.3f %7s %5.0f%%  %s"
                  % (name, truth_cm, mu * 100, err * 100,
                     ("%+.1f%%" % (err / (truth_cm / 100) * 100)) if truth_cm else "—",
                     100.0 * rej / a.seeds,
                     desc + ("" if len(ests) < 2 else "   σ=%.3f" % (sd * 100))))
        else:
            print("%-15s %6.2f %9s %8s %7s %5.0f%%  %s"
                  % (name, truth_cm, "—", "—", "—", 100.0, desc))


if __name__ == "__main__":
    main()
