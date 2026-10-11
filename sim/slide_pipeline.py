"""真机滑动测量分析：按钮边界 + 算法审计

这是 sim/pipeline.py 的「真机版」：不再有真值，边界由 App 的触摸事件给定。
它做四件事：

  1. 时间对齐：把 meta.txt 里 DOWN/UP 的 recv_elapsed_ns 落到 acc.csv 的时间轴上
  2. 校准：用 DOWN 之前的静止段估陀螺零偏 / 重力参考 / 噪声地板
  3. 审计：检查 DOWN/UP 时刻附近手机是否真的静止（按钮方案唯一的失效模式）
  4. 出距离：按钮窗口版完整管线，并与「检测器自动定位」对照

用法：
    python slide_pipeline.py <采集目录>
不给参数则取 ../android/captures/slideruler 下最新的一个。
"""

import csv
import math
import os
import re
import statistics
import sys

import pipeline as P
from pipeline import G, DeviceProfile, Params, vdot, vnorm, vscale, vsub, vunit

CAPTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "android", "captures", "slideruler")

AUDIT_HALF_S = 0.3      # 审计窗半径
CAL_MAX_S = 2.0         # 校准窗上限


# --------------------------------------------------------------------------
def read_pairs(path):
    ts, vals = [], []
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) < 4:
                continue
            ts.append(int(row[0]))
            vals.append((float(row[1]), float(row[2]), float(row[3])))
    return ts, vals


def read_marks(path):
    marks = {}
    rx = re.compile(r"(DOWN|UP)\s+event=(\d+) ms.*recv_elapsed=(\d+) ns.*输入延迟=(-?\d+) ms")
    with open(path, encoding="utf-8") as f:
        for line in f:
            m = rx.search(line)
            if m:
                marks[m.group(1)] = {"event_ms": int(m.group(2)),
                                     "recv_ns": int(m.group(3)),
                                     "latency_ms": int(m.group(4))}
    return marks


def locate(ts, target_ns):
    """二分找到第一个 >= target 的样本下标"""
    lo, hi = 0, len(ts) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if ts[mid] < target_ns:
            lo = mid + 1
        else:
            hi = mid
    return lo


# --------------------------------------------------------------------------
def analyse(d):
    print("=" * 92)
    print(f"采集目录  {os.path.basename(d)}")
    print("=" * 92)

    acc_ts, acc = read_pairs(os.path.join(d, "acc.csv"))
    gyr_ts, gyr = read_pairs(os.path.join(d, "gyr.csv"))
    marks = read_marks(os.path.join(d, "meta.txt"))

    n = min(len(acc), len(gyr))
    acc, gyr = acc[:n], gyr[:n]

    dts = [acc_ts[i + 1] - acc_ts[i] for i in range(min(len(acc_ts), n) - 1)]
    dt = statistics.median(dts) / 1e9
    rate = 1.0 / dt
    dur = n * dt
    print(f"样本 acc {len(acc)} / gyr {len(gyr)}   实测 dt {dt*1000:.4f} ms ({rate:.2f} Hz)   "
          f"时长 {dur:.2f} s")
    if not marks:
        print("meta.txt 里没有触摸事件 —— 这不是一次「按住式」采集。")
        return

    # ---- 1. 时间对齐 ----
    t0 = acc_ts[0]
    down_ns = marks["DOWN"]["recv_ns"]
    up_ns = marks.get("UP", {}).get("recv_ns")
    lo = locate(acc_ts, down_ns)
    hi = locate(acc_ts, up_ns) if up_ns else None
    print()
    print("1 · 时间对齐")
    print(f"  DOWN  相对首样本 {(down_ns-t0)/1e6:>10.2f} ms  → 样本 #{lo}"
          f"   事件→送达 {marks['DOWN']['latency_ms']} ms")
    if hi is not None:
        print(f"  UP    相对首样本 {(up_ns-t0)/1e6:>10.2f} ms  → 样本 #{hi}"
              f"   事件→送达 {marks['UP']['latency_ms']} ms")
        print(f"  按住跨度 {(up_ns-down_ns)/1e6:.2f} ms = {(hi-lo)*dt*1000:.2f} ms 的采样")
        # 事件时刻与最近采样点的错位（按钮方案里，这一点错位被静止段吸收）
        off_lo = (down_ns - acc_ts[max(0, lo - 1)]) / 1e6
        print(f"  DOWN 落在两个采样点之间：距上一采样 {off_lo:.2f} ms"
              f"（起点边界最坏差一个采样点 = {dt*1000:.2f} ms）")

    # ---- 2. 校准（DOWN 之前的静止段）----
    cal_n = max(20, min(int(CAL_MAX_S / dt), lo))
    cal_lo = lo - cal_n
    b_g = tuple(statistics.fmean(gyr[i][a] for i in range(cal_lo, lo)) for a in range(3))
    a_mean = tuple(statistics.fmean(acc[i][a] for i in range(cal_lo, lo)) for a in range(3))
    a_sig = max(statistics.pstdev([acc[i][a] for i in range(cal_lo, lo)]) for a in range(3))
    g_sig = max(statistics.pstdev([gyr[i][a] for i in range(cal_lo, lo)]) for a in range(3))
    prof = DeviceProfile(dt, b_g, vunit(a_mean), vnorm(a_mean), a_sig, g_sig, cal_n)
    print()
    print("2 · 校准（用 DOWN 之前 %.2f s 的静止段）" % (cal_n * dt))
    print(f"  {prof.summary()}")

    # ---- 3. 姿态 + 去重力 ----
    p = Params(static_mode="adaptive", assumed_dt=dt)
    truth = {"N": n - 1, "acc": acc, "gyr": gyr, "dt": dt,
             "n_pre": lo, "n_mot": (hi - lo) if hi else 0}
    qs, _ = P.solve_attitude(truth, prof.gyro_bias, prof.gravity_ref, p, dt=dt)
    # g_mag 必须用「本机静止段实测的重力模长」，不能写死 9.81。
    # K40 实测 ‖g‖=9.911，与 9.81 差 0.10 m/s² —— 正好等于静止段残余 |a_lin| 的量级，
    # 会把 PCA 主轴从"运动方向"带偏到"重力方向"（慢速滑动时会直接输出 ~0）。
    a_lin = P.remove_gravity(truth, qs, g_mag=prof.gravity_mag)

    # ---- 4. 审计 ----
    w = max(1, int(AUDIT_HALF_S / dt))

    def seg_max(idx, stats_lo=None, stats_hi=None):
        a_ = max(0, idx - w)
        b_ = min(len(a_lin), idx + w + 1)
        return max(vnorm(a_lin[i]) for i in range(a_, b_))

    base_lo, base_hi = max(0, cal_lo), lo
    base = max(vnorm(a_lin[i]) for i in range(base_lo, base_hi))
    thr = 1.6 * base + 0.05

    print()
    print("3 · 审计（按钮方案唯一的失效模式是「按下时手机已在动 / 松手时还没停」）")
    print(f"  静止段基准 max‖a_lin‖ = {base:.4f} m/s²   自适应阈值 = {thr:.4f} m/s²")
    ok_down = seg_max(lo) < thr
    print(f"  DOWN ±{AUDIT_HALF_S}s   max‖a_lin‖ = {seg_max(lo):>8.4f}   "
          f"{'✓ 静止' if ok_down else '✗ 可疑（按下时手机似乎在动）'}")
    if hi is not None:
        ok_up = seg_max(hi) < thr
        print(f"  UP   ±{AUDIT_HALF_S}s   max‖a_lin‖ = {seg_max(hi):>8.4f}   "
              f"{'✓ 静止' if ok_up else '✗ 可疑（松手时手机似乎还没停）'}")

    # ---- 5. 距离 ----
    # 静止段残余零偏量级 b（投影到主轴后再算，见下）与噪声地板 σ_a：
    cal_range = range(max(0, lo - int(CAL_MAX_S / dt)), lo)
    sigma_a = prof.acc_sigma   # 逐轴噪声地板取最大（不能用三轴混在一起算：z 轴含 9.8）
    # 残余零偏量级：校准段 a_lin 的均值模长（用实测重力模长去重力后，这个值应远小于 0.02）
    b_vec = tuple(statistics.fmean(a_lin[i][a] for i in cal_range) for a in range(3))
    bias_bound = vnorm(b_vec)
    print()
    print(f"  静止段残余零偏 |b| = {bias_bound:.4f} m/s²   acc 噪声 σ = {sigma_a:.4f} m/s² "
          f"（b·T 决定条件 ZUPT 的阈值）")

    def distance(lo_i, hi_i):
        """返回 (距离, 裸积分, v_end, 阈值, 是否去趋势, λ1/λ2)"""
        u, lam1, lam2 = P.pca_axis(a_lin[lo_i:hi_i + 1])
        axis = u if vdot(u, (1.0, 0.0, 0.0)) >= 0 else vscale(u, -1.0)
        a1 = [vdot(a_lin[i], axis) for i in range(len(a_lin))]
        _, p_raw = P.integrate(a1, lo_i, hi_i, dt=dt)
        # 残余零偏投影到主轴（1D）
        b1 = statistics.fmean(a1[i] for i in cal_range)
        dist, v_end, used, thr = P.integrate_cond_zupt(
            a1, lo_i, hi_i, dt, bias_bound=abs(b1), sigma=sigma_a)
        return dist, abs(p_raw[-1]), v_end, thr, used, (lam1 / lam2 if lam2 > 0 else float("inf"))

    print()
    print("5 · 距离（条件 ZUPT：只有 |v_end| 大到零偏解释不了时才去趋势）")
    if hi is not None and hi > lo + 5:
        dist, raw, v_end, thr, used, lr = distance(lo, hi)
        T = (hi - lo) * dt
        print(f"  按钮窗口  [{lo}, {hi}]  时长 {T:.2f} s")
        print(f"    λ1/λ2 = {lr:.1f}    终点速度 v_end = {v_end:+.4f} m/s    阈值 = {thr:.4f} m/s"
              f"    → {'去趋势' if used else '裸积分'}")
        print(f"    ★ 距离 = {dist*100:.2f} cm      （裸积分 {raw*100:.2f} cm）")
        if used:
            print(f"      去趋势扣掉 {abs(raw-dist)*100:.2f} cm（小 v_end 时≈|v_end|·T/2={abs(v_end)*T/2*100:.2f} cm）")
        # 把窗口两端各外扩 200ms，验证「余量无关性」
        m = int(0.2 / dt)
        lo2, hi2 = max(0, lo - m), min(len(a_lin) - 1, hi + m)
        d2, raw2, v2, thr2, u2, _ = distance(lo2, hi2)
        print(f"    外扩 200ms [{lo2}, {hi2}]  {d2*100:.2f} cm   "
              f"→ 与上者差 {abs(d2-dist)*1000:.2f} mm（应 ≈0）")
    else:
        print("  没有可用的 UP 边界")

    # ---- 对照：检测器自动定位 ----
    flags = P.detect_static(truth, b_g, p, a_lin=a_lin, profile=prof)
    dlo, dhi = P.motion_window(flags, truth, p)
    print()
    print("6 · 对照：静止检测器自动定位的窗口")
    if dlo is None or dhi is None:
        print("  检测器没找到运动段（本次手机确实没动，符合预期）")
    else:
        dd, draw, dv, dthr, dused, _ = distance(dlo, dhi)
        print(f"  [{dlo}, {dhi}]  {dd*100:.2f} cm（裸积分 {draw*100:.2f}）   "
              f"v_end={dv:+.4f} {'去趋势' if dused else '裸积分'}   "
              f"（按钮窗 {lo}~{hi}，检测器窗 {dlo}~{dhi}）")
    print()


def main():
    if len(sys.argv) > 1:
        dirs = [sys.argv[1]]
    else:
        if not os.path.isdir(CAPTURES):
            print(f"找不到采集目录 {CAPTURES}")
            return 1
        all_d = sorted([os.path.join(CAPTURES, x) for x in os.listdir(CAPTURES)
                        if os.path.isdir(os.path.join(CAPTURES, x))])
        dirs = all_d[-1:]
    for d in dirs:
        analyse(d)
    return 0


if __name__ == "__main__":
    sys.exit(main())
