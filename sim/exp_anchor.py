#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""实验：锚点到底该用多长？

背景（ChatGPT 代码审查第 6 条引出来的）
------------------------------------------------------------------
`STATIC_MIN_S = 0.20` 声明了「静止段至少 0.20 秒」，但两端都没落实。
更要命的是：窗口 = 运动段 ± MARGIN_S(0.15)，而静止段被裁进窗口之后
**剩下的恰好就是那段余量** —— 实测有效锚点恒为 0.143 s（15 个样本），
一次都没多过。所以：

  · 照字面执行 STATIC_MIN_S = 0.20，7 条采集会**全部被拒**（0.143 < 0.20）；
  · 「双锚点把静止段全部样本拉平，比只钉两端稳 √N」这句话里的 N 恒等于 15，
    √N 的收益被余量卡死了 —— 设计文档里那句话是错的。

要回答的物理问题
------------------------------------------------------------------
静止段延长窗口到底安不安全？

反对的理由（我原先写在代码注释里的）：静止段对真实位移零贡献，
但对误差是完整加成，误差 ∝ T²。
支持的理由：那个 ∝T² 是**裸积分**的账。ZUPT 之后，静止段里速度按
≈b·t 线性增长，恰恰是拟合 A、B 最干净的数据；延长窗口等于白送样本。

两边都有道理，所以只能看数据。本脚本在三种窗口上算同一个数：
  margin   现设计：运动段 ± MARGIN_S（锚点 15 样）
  motion   贴着运动段：不放余量（锚点 ≈0，退化，仅作参照）
  full     外扩到完整静止段（锚点 = 整段静止）

并同时报告锚点内速度的散布 —— 那是"锚得牢不牢"的直接证据。

用法：python exp_anchor.py
"""
import os
import re
import statistics

import measure2 as M

CAP_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "android", "captures", "slideruler")

# 用户肉眼估的真值，本身 ±0.5 cm（±3.3%）
TRUTH = {
    "20261011_100444": 15.0,
    "20261011_100459": 15.0,
    "20261011_100529": 15.0,
}


def load(d):
    acc_ts, acc = M.read_pairs(os.path.join(d, "acc.csv"))
    gyr_ts, gyr = M.read_pairs(os.path.join(d, "gyr.csv"))
    meta = open(os.path.join(d, "meta.txt"), encoding="utf-8").read()
    dn = int(re.search(
        r"DOWN\s+event=\d+ ms\s+recv_uptime=\d+ ms\s+recv_elapsed=(\d+) ns", meta).group(1))
    up = int(re.search(
        r"UP\s+event=\d+ ms\s+recv_uptime=\d+ ms\s+recv_elapsed=(\d+) ns", meta).group(1))
    n = min(len(acc), len(gyr))
    acc, gyr = acc[:n], gyr[:n]
    dt = statistics.median([acc_ts[i + 1] - acc_ts[i] for i in range(n - 1)]) / 1e9

    def locate(t):
        lo, hi = 0, n - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if acc_ts[mid] < t:
                lo = mid + 1
            else:
                hi = mid
        return lo

    return acc, gyr, locate(dn), locate(up), dt


def run_containing(flags, i, n):
    a = i
    while a > 0 and flags[a - 1]:
        a -= 1
    b = i + 1
    while b < n and flags[b]:
        b += 1
    return (a, b)


def zupt_full(a_lin, w0, w1, flags, n, dt):
    """用完整静止段当锚点，返回 (距离, 左锚点样数, 右锚点样数, 段内速度散布)"""
    dist, vc, p, drift, spread = M.zupt_two_anchor(
        a_lin, w0, w1, run_containing(flags, w0, n), run_containing(flags, w1, n), dt)
    return dist, vc, p, drift, spread


def anchor_stats(vc, w0, run, dt):
    """锚点（裁进窗口后）内的速度：样本数与模长标准差 —— 散布小=锚得牢"""
    a = max(run[0], w0)
    b = min(run[1], w0 + len(vc))
    if b - a < 1:
        return 0, float("nan")
    mags = [M.vnorm(vc[i - w0]) for i in range(a, b)]
    return b - a, (statistics.pstdev(mags) if len(mags) > 1 else 0.0)


def main():
    print("窗口三选一：锚点越长，段内速度均值越稳（√N），但窗口也越长（∝T²）")
    print()
    print("%-16s %6s | %-24s | %-24s | %-24s" % (
        "采集", "真值", "margin (现设计)", "motion (无余量)", "full (外扩到整段静止)"))
    print("-" * 118)

    rows = []
    for nm in sorted(os.listdir(CAP_ROOT)):
        d = os.path.join(CAP_ROOT, nm)
        if not os.path.isdir(d):
            continue
        try:
            acc, gyr, dn, up, dt = load(d)
        except Exception:
            continue
        r = M.measure(acc, gyr, dn, up, dt=dt)
        if not (r and r["ok"] and not r.get("noMotion")):
            continue

        flags = r["flags"]
        m0, m1 = r["motion"]
        run0, run1 = r["run0"], r["run1"]
        n = len(flags)

        # 重算 a_lin（measure 已经把 gb2/g_ref/gMag 吐出来了，直接用）
        qs = M.solve_attitude_v2(acc, gyr, r["gyroBias"], r["gRef"], flags, dt)
        a_lin = M.to_linear(acc, qs, r["gMag"])

        marg = max(1, int(M.MARGIN_S / dt))
        variants = [
            ("margin", max(0, m0 - marg), min(n - 1, m1 + marg)),
            ("motion", m0, m1),
            ("full", run0[0], run1[1] - 1),
        ]

        cells = []
        for label, w0, w1 in variants:
            if w1 - w0 < 5:
                cells.append((float("nan"), 0, 0, float("nan")))
                continue
            dist, vc, p, drift, spread = zupt_full(a_lin, w0, w1, flags, n, dt)
            # 锚点（裁进窗口之后）的样本数
            a0 = max(run_containing(flags, w0, n)[0], w0)
            b0 = min(run_containing(flags, w0, n)[1], w1 + 1)
            a1 = max(run_containing(flags, w1, n)[0], w0)
            b1 = min(run_containing(flags, w1, n)[1], w1 + 1)
            nmin = min(b0 - a0, b1 - a1)
            mags = [M.vnorm(vc[i - w0]) for i in range(a0, min(b0, w0 + len(vc)))]
            sd = statistics.pstdev(mags) if len(mags) > 1 else float("nan")
            cells.append((dist, nmin, (w1 - w0), sd))

        t = TRUTH.get(nm)
        err = ""
        if t:
            err = "  (full %.2f%% / margin %.2f%%)" % (
                (cells[2][0] - t / 100) / (t / 100) * 100,
                (cells[0][0] - t / 100) / (t / 100) * 100)
        print("%-16s %6s | %s%s" % (
            nm, t if t else "-",
            "  ".join("%6.2fcm 锚%3d样 T%.2fs" % (c[0] * 100, c[1], c[2] * dt)
                      for c in cells),
            err))
        rows.append((nm, t, cells))

    print()
    print("怎么读")
    print("  锚 N 样 = 左右两侧较短那个锚点、裁进窗口后剩下的样本数。")
    print("  margin 那一列恒为 15 样 —— 这就是「余量卡死 √N」的证据。")
    print("  full 那一列是外扩到完整静止段后的锚点数；若它的距离与 margin 接近，")
    print("  说明静止段延长窗口是安全的（∝T² 那一派错了）；若差得远，说明有别的账。")
    print()
    print("有真值的三条，full 与 margin 的偏差：")
    for nm, t, cells in rows:
        if not t:
            continue
        dm = (cells[0][0] - t / 100) / (t / 100) * 100
        df = (cells[2][0] - t / 100) / (t / 100) * 100
        print("  %s  真值 %.0f cm   margin %+.2f%% (%d 样)   full %+.2f%% (%d 样)   "
              "→ %s" % (nm, t, dm, cells[0][1], df, cells[2][1],
                        "full 更好" if abs(df) < abs(dm) else "margin 更好"))


if __name__ == "__main__":
    main()
