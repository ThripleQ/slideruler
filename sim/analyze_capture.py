#!/usr/bin/env python3
"""
真机采集体检 —— android/ 落下的 acc.csv / gyr.csv 的原始时间序列二次分析。

体检 A 的汇总报告（meta.txt）回答不了这三类问题，而它们恰好决定项目精度：

  1. 时间戳节奏
     报告只给"间隔中位 / P95 / max"。但如果 dt 恒定不变，那说明时间戳是**框架合成的**
     （传感器批次化 batch → SensorEventQueue 均匀铺开），而不是硬件真值。这很关键：
     合成时间戳意味着每个样本的真实时刻有一个"批次化级别"的定时不确定度，而仿真已经
     证明**起点边界定位是全项目的主误差项**（起点晚 10 ms → -1.65%）。

  2. 噪声 / 漂移
     整体 σ 会被起步瞬态（手指点按钮那一下）和慢漂污染。这里给 分窗 σ / 去漂 σ /
     一阶差分 σ / Allan 偏差，用来把"白噪声"和"零偏不稳定性"分开。

  3. 零偏与标度
     三轴均值、模长偏离 9.81、隐含倾角；两次独立重复对照看重复性。

用法：
    python analyze_capture.py                       # 自动扫 ../android/captures/slideruler/*
    python analyze_capture.py <capture_dir> [...]
"""
import csv
import math
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ROOT = os.path.normpath(os.path.join(HERE, "..", "android", "captures", "slideruler"))

G_NOMINAL = 9.81
AXES = ("x", "y", "z")


# ------------------------------------------------------------------ 基础统计

def mean(x):
    return sum(x) / len(x) if x else 0.0


def std(x):
    if len(x) < 2:
        return 0.0
    m = mean(x)
    return math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1))


def pct(x, p):
    if not x:
        return 0.0
    s = sorted(x)
    i = min(len(s) - 1, max(0, int(round(p / 100.0 * (len(s) - 1)))))
    return s[i]


def detrend(x):
    """最小二乘去线性趋势（去掉慢漂，剩下的更接近传感器白噪声）。"""
    n = len(x)
    if n < 3:
        return list(x)
    mx = (n - 1) / 2.0
    my = mean(x)
    num = sum((i - mx) * (x[i] - my) for i in range(n))
    den = sum((i - mx) ** 2 for i in range(n))
    b = num / den if den else 0.0
    a = my - b * mx
    return [x[i] - (a + b * i) for i in range(n)]


def diff_std(x):
    """一阶差分 / sqrt(2)：只由样本间高频起伏决定，对任何慢漂免疫。"""
    if len(x) < 3:
        return 0.0
    d = [x[i + 1] - x[i] for i in range(len(x) - 1)]
    return std(d) / math.sqrt(2.0)


def allan(x, dt, tau):
    """重叠 Allan 偏差（非重叠块版本）。

    分辨噪声类型的标准工具：
      - ADEV ∝ τ^-0.5  → 白噪声（角随机游走/RRW 的前段）
      - ADEV 在某个 τ 处出现平台 → 零偏不稳定性 bias instability（sim 里的 gyro_instab）
      - ADEV ∝ τ^+0.5 → 速率随机游走
    """
    m = max(1, int(round(tau / dt)))
    nb = len(x) // m
    if nb < 3:
        return None
    means = [mean(x[i * m:(i + 1) * m]) for i in range(nb)]
    d = [means[i + 1] - means[i] for i in range(nb - 1)]
    return math.sqrt(0.5 * mean([v * v for v in d]))


def lag1_autocorr(x):
    if len(x) < 3:
        return 0.0
    m = mean(x)
    num = sum((x[i] - m) * (x[i + 1] - m) for i in range(len(x) - 1))
    den = sum((v - m) ** 2 for v in x)
    return num / den if den else 0.0


# ------------------------------------------------------------------ 读取

def load(path):
    ts = []
    ax = ([], [], [])
    with open(path, newline="") as f:
        rd = csv.reader(f)
        next(rd, None)
        for row in rd:
            if len(row) < 4:
                continue
            try:
                t = int(row[0])
                v = (float(row[1]), float(row[2]), float(row[3]))
            except ValueError:
                continue
            ts.append(t)
            for k in range(3):
                ax[k].append(v[k])
    return ts, ax


# ------------------------------------------------------------------ 三个分析块

def report_clock(name, ts):
    print(f"  ── 时间戳节奏（{name}）──")
    n = len(ts)
    if n < 3:
        print("    样本太少")
        return
    span = (ts[-1] - ts[0]) / 1e9
    d = [ts[i + 1] - ts[i] for i in range(n - 1)]
    cnt = Counter(d)
    dm = mean(d) / 1e6
    print(f"    样本 {n}   时长 {span:.3f} s   平均间隔 {dm:.6f} ms   ({1000.0 / dm:.3f} Hz)")
    print(f"    不同的 dt 取值个数：{len(cnt)}   （=1 即完全恒定 → 框架合成时间戳）")
    for v, c in cnt.most_common(5):
        print(f"      {v / 1e6:>14.6f} ms  ×{c}   ({c / len(d) * 100:.1f}%)")
    lo, hi = min(d), max(d)
    print(f"    dt 范围 {(hi - lo) / 1e6:.6f} ms   抖动(σ) {std(d) / 1e6:.6f} ms"
          f"   P95 {(pct(d, 95)) / 1e6:.4f} ms   max {hi / 1e6:.4f} ms")
    # 批次化指纹：是否存在若干个"正常间隔"的整数倍
    base = cnt.most_common(1)[0][0]
    multi = sorted({round(v / base) for v in cnt if v > base * 1.5})
    if multi:
        print(f"    ⚠ 出现基本间隔的整数倍：{multi} → 存在批次化/丢样拼接")
    else:
        print(f"    ✓ 没有基本间隔的整数倍 → 不存在批次化拼接")
    # 时间戳是否落在固定网格上
    frac = {(ts[i] - ts[0]) % base for i in range(0, n, max(1, n // 50))}
    print(f"    相对首样本取模基本间隔的余数种类：{len(frac)}"
          f"  （1~2 种 = 完美网格 / 真实抖动）")


def report_noise(name, ax, dt):
    print(f"  ── 噪声 / 漂移（{name}）──")
    print(f"    {'轴':<4}{'均值':>11}{'整体σ':>11}{'去漂σ':>11}{'一阶差分σ':>12}{'lag1自相关':>12}")
    for k in range(3):
        x = ax[k]
        print(f"    {AXES[k]:<4}{mean(x):>11.5f}{std(x):>11.5f}"
              f"{std(detrend(x)):>11.5f}{diff_std(x):>12.5f}{lag1_autocorr(x):>12.4f}")
    n = len(ax[0])
    step = max(1, n // 6)
    print(f"    分窗均值（每窗 ≈{step * dt:.1f} s，看是否慢漂）")
    hdr = "      " + "".join(f"{f'w{i}':>11}" for i in range(min(6, n // step + 1)))
    print(hdr)
    for k in range(3):
        vals = []
        for i in range(0, n, step):
            seg = ax[k][i:i + step]
            if len(seg) < 5:
                continue
            vals.append(mean(seg))
        print(f"      {AXES[k]:<4}" + "".join(f"{v:>11.5f}" for v in vals[:6]))
    print(f"    Allan 偏差（分辨白噪声 τ^-0.5 vs 零偏不稳定平台）")
    print(f"      {'τ(s)':>7}" + "".join(f"{a:>12}" for a in AXES))
    for tau in (0.1, 0.5, 1.0, 2.0, 5.0, 10.0):
        row = []
        for k in range(3):
            a = allan(detrend(ax[k]), dt, tau)
            row.append("--" if a is None else f"{a:.6f}")
        print(f"      {tau:>7.1f}" + "".join(f"{v:>12}" for v in row))


def report_bias(name, ax, unit_note=""):
    v = [mean(ax[k]) for k in range(3)]
    mag = math.sqrt(sum(t * t for t in v))
    print(f"  ── 零偏 / 模长（{name}）──{unit_note}")
    print(f"    均值向量 ({v[0]:+.5f}, {v[1]:+.5f}, {v[2]:+.5f})   模长 {mag:.5f}")
    if mag > 1.0:  # 加速度计
        dev = mag - G_NOMINAL
        tilt = math.degrees(math.acos(min(1.0, abs(v[2]) / mag)))
        hx, hy = v[0], v[1]
        print(f"    模长偏离 9.81 = {dev:+.5f} m/s²  ({dev / G_NOMINAL * 100:+.3f}%)")
        print(f"    z 轴与重力夹角 {tilt:.3f}°   水平分量 ({hx:+.5f}, {hy:+.5f})"
              f"  水平倾角 {math.degrees(math.atan2(math.hypot(hx, hy), abs(v[2]))):.3f}°")
    else:  # 陀螺
        print(f"    模长 {mag:.6f} rad/s = {math.degrees(mag):.4f} °/s")


def analyse(d, out=None):
    print("=" * 78)
    print(f"采集目录  {os.path.basename(d)}")
    meta = os.path.join(d, "meta.txt")
    if os.path.exists(meta):
        with open(meta, encoding="utf-8") as f:
            for line in f:
                if line.startswith(("时间", "设备", "样本数")) or "样本 / gyr" in line:
                    print("   " + line.rstrip())
    acc_ts, acc = load(os.path.join(d, "acc.csv"))
    gyr_ts, gyr = load(os.path.join(d, "gyr.csv"))
    if not acc_ts:
        print("  空")
        return
    dt = (acc_ts[-1] - acc_ts[0]) / 1e9 / (len(acc_ts) - 1)
    report_clock("加速度计", acc_ts)
    report_clock("陀螺仪", gyr_ts)
    print()
    report_bias("加速度计", acc, "  [m/s²]")
    report_bias("陀螺仪", gyr, "  [rad/s]")
    print()
    report_noise("加速度计", acc, dt)
    print()
    report_noise("陀螺仪", gyr, dt)
    print()
    if out is not None:
        n = len(acc[0])
        span = (acc_ts[-1] - acc_ts[0]) / 1e9
        out.append({
            "dir": os.path.basename(d),
            "acc_mean": tuple(mean(acc[k]) for k in range(3)),
            "acc_mag": math.sqrt(sum(mean(acc[k]) ** 2 for k in range(3))),
            "gyr_mag": math.sqrt(sum(mean(gyr[k]) ** 2 for k in range(3))),
            "acc_sigma": tuple(std(acc[k]) for k in range(3)),
            "rate": n / span if span > 0 else 0.0,
        })


def report_cross(items):
    """跨采集对照：判定"模长随姿态变化"还是"等比标度误差"。

    物理依据：静止时真实比力向量的模长恒等于当地 g，**与姿态无关**（姿态只旋转向量，
    不改长度）。所以：
      - 模长在各姿态下基本不变 → 误差更像**等比标度**（三轴同比例放大）
      - 模长随姿态明显变化   → 必然是**零偏/非等比**（有个与体轴绑定的常向量）
    """
    print("=" * 78)
    print("跨采集对照 —— 姿态相关性检验")
    print(f"  {'目录':<18}{'acc 均值 (x, y, z)':<38}{'模长':>9}{'偏离':>10}{'gyr模长':>10}")
    for it in items:
        m = it["acc_mean"]
        print(f"  {it['dir']:<18}"
              f"({m[0]:+.4f}, {m[1]:+.4f}, {m[2]:+.4f})".ljust(38)
              + f"{it['acc_mag']:>9.4f}{it['acc_mag'] - G_NOMINAL:>+10.4f}"
                f"{it['gyr_mag']:>10.5f}")
    mags = [it["acc_mag"] for it in items]
    spread = max(mags) - min(mags)
    print(f"\n  模长极差 = {spread:.4f} m/s²  （姿态只旋转向量、不改长度，理论上应≈0）")
    if spread > 0.02:
        print("  ⇒ 模长随姿态变化：**排除等比标度误差**，必然存在与体轴绑定的零偏。")
    else:
        print("  ⇒ 模长与姿态无关：更像等比标度误差（三轴同比例）。")

    groups = {}
    for it in items:
        m = it["acc_mean"]
        k = max(range(3), key=lambda i: abs(m[i]))
        groups.setdefault(k, {})[1 if m[k] >= 0 else -1] = it
    print("\n  按主导轴配对（同轴有正反两姿态时可粗解该轴零偏/标度）")
    for k in sorted(groups):
        g = groups[k]
        if 1 in g and -1 in g:
            mp, mn = g[1]["acc_mag"], g[-1]["acc_mag"]
            b = (mp - mn) / 2.0
            kg = (mp + mn) / 2.0
            print(f"    {AXES[k]} 轴  |f|₊={mp:.4f}  |f|₋={mn:.4f}"
                  f"   ⇒ 沿轴零偏 ≈ {b:+.4f} m/s²   等效 k·g ≈ {kg:.4f}"
                  f"（{kg / G_NOMINAL - 1:+.3%}）")
        else:
            print(f"    {AXES[k]} 轴  只有{'正' if 1 in g else '负'}向姿态，无法分离")
    print("\n  注：上面两式忽略了倾斜与横向分量的耦合，是粗估；严格解需要受控多姿态标定"
          "（每个轴正反各一次），或直接做已知距离的端到端标定。")
    print()


def main():
    args = sys.argv[1:]
    if args:
        dirs = args
    else:
        if not os.path.isdir(DEFAULT_ROOT):
            print(f"找不到默认目录 {DEFAULT_ROOT}")
            return 1
        dirs = sorted(
            os.path.join(DEFAULT_ROOT, x)
            for x in os.listdir(DEFAULT_ROOT)
            if os.path.isdir(os.path.join(DEFAULT_ROOT, x))
        )
    if not dirs:
        print("没有采集目录")
        return 1
    collected = []
    for d in dirs:
        analyse(d, collected)
    if len(collected) > 1:
        report_cross(collected)
    return 0


if __name__ == "__main__":
    sys.exit(main())
