#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
实验 14：测量区间到底该由谁决定 —— 静止检测，还是按钮？

背景（用户 2026-10-11 的原话）
----------------------------
「滑动过程中会被判定静止，导致长度被拆，最终给出的是一部分长度。」
「为什么一定要判定静止状态？以我的物理知识，只要从静止启动就好了。
  只要静止启动，就认为按钮按住时采集的是一段完整的位移，按钮松开后那不是测量结束了吗？
  及时松开时还在移动，但数据不是已经到了吗，按钮松开后为什么还要判定静止呢？」

两条其实是同一个问题：**静止检测被同时用来做两件不同的事**——
  (1) 定义「测量区间」（哪一段要积分）   ← 用户反对的就是这个用法
  (2) 提供「误差定标锚点」（残余误差是多少） ← 这才是它该干的活
(1) 一旦用了静止检测，匀速滑动就会被判成静止（加速度计分不出匀速和静止，
    这是伽利略相对性原理，不是算法的锅），于是一次完整滑动被切成好几段，
    只取其中一段 → 距离被腰斩。真机证据：122416 一次滑动被切成 6 段，
    只取了 0.29 s → 报 0.59 cm；122407 切成 5 段，只取 0.76 s → 报 3.00 cm。

所以本实验把这两件事拆开，比三种方案：
  A 现役   ：区间 = 静止检测切出的运动段 ± 余量，双锚点
  B 按钮窗 ：区间 = [DOWN, UP] ± 余量（**完全不做运动分段**），双锚点
  C 只起点 ：区间 = [DOWN, UP] ± 余量，只用按下前那段静止估一个常值残差减掉
             —— 这就是用户提的方案：「静止启动就够了」
C 相对 B 差多少，就是「终点静止到底值多少钱」的实测答案。

用法：python exp_interval.py
"""

import os
import re
import statistics

import measure2 as M

CAP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "android", "captures", "slideruler")

# 用户肉眼估的真值（cm）。本身 ±0.5 cm ≈ ±3.3%。
TRUTH = {
    "20261011_100444": 15.0,
    "20261011_100459": 15.0,
    "20261011_100529": 15.0,
}


def load(d):
    at, acc = M.read_pairs(os.path.join(d, "acc.csv"))
    gt, gyr = M.read_pairs(os.path.join(d, "gyr.csv"))
    n = min(len(acc), len(gyr))
    acc, gyr = acc[:n], gyr[:n]
    dt = statistics.median([at[i + 1] - at[i] for i in range(n - 1)]) / 1e9
    meta = open(os.path.join(d, "meta.txt"), encoding="utf-8").read()
    m = re.search(r"DOWN\s+event=\d+ ms\s+recv_uptime=\d+ ms\s+recv_elapsed=(\d+) ns", meta)
    if not m:
        return None
    dn_ns = int(m.group(1))
    mu = re.search(r"UP\s+event=\d+ ms\s+recv_uptime=\d+ ms\s+recv_elapsed=(\d+) ns", meta)
    if not mu:
        return None
    up_ns = int(mu.group(1))

    def locate(t):
        lo, hi = 0, n - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if at[mid] < t:
                lo = mid + 1
            else:
                hi = mid
        return lo

    return acc, gyr, dt, locate(dn_ns), locate(up_ns), n


def run_containing(flags, i, n):
    a = i
    while a > 0 and flags[a - 1]:
        a -= 1
    b = i + 1
    while b < n and flags[b]:
        b += 1
    return (a, b)


def variants(d):
    """返回 A / B / C 三个方案的距离（米）与说明。"""
    got = load(d)
    if got is None:
        return None
    acc, gyr, dt, dn, up, n = got
    if up <= dn or n < 60:
        return None

    m = max(1, int(M.MARGIN_S / dt))
    # ---- 静止检测（与现役完全相同的判据）----
    pre = max(20, min(n // 4, int(0.5 / dt)))
    sd0 = M.rolling_max_std(acc, dt)
    iq = min(range(n), key=lambda i: sd0[i])
    qa, qb = max(0, iq - pre // 2), min(n, iq + pre // 2 + 1)
    gb = tuple(statistics.fmean(gyr[i][a] for i in range(qa, qb)) for a in range(3))
    flags = M.static_flags(acc, gyr, gb, dt)[0]
    min_run = max(3, int(M.STATIC_MIN_S / dt))

    # ---- A：现役 ----
    ra = M.measure(acc, gyr, dn, up, dt=dt)
    da = ra["distanceM"] if (ra["ok"] and not ra.get("noMotion")) else None

    # ---- B / C：窗口 = [DOWN, UP] ± 余量，不做任何运动分段 ----
    w0 = max(0, dn - m)
    w1 = min(n - 1, up + m)
    nb = {}
    run0 = run_containing(flags, w0, n)
    run1 = run_containing(flags, w1, n)
    if not flags[w0]:
        nb["why"] = "按下那一刻手机不在静止状态"
    elif not flags[w1]:
        nb["why"] = "松手那一刻手机不在静止状态"
    elif run0[0] == run1[0]:
        nb["why"] = "按下与松手落在同一段静止里（没动）"
    elif min(run0[1] - run0[0], run1[1] - run1[0]) < min_run:
        nb["why"] = "两端的静止段有一侧太短"

    db = dc = None
    cal = None
    if not nb:
        cal = (max(run0[0], dn - int(M.CAL_MAX_S / dt)), dn)
        if (cal[1] - cal[0]) * dt < M.CAL_MIN_S:
            nb["why"] = "按下前的静止段太短（%.2f s）" % ((cal[1] - cal[0]) * dt)
        else:
            ca, cb = cal
            fr = sum(1 for i in range(ca, cb) if flags[i]) / (cb - ca)
            if fr < M.CAL_STATIC_FRAC:
                nb["why"] = "校准段不够静止（%.0f%%）" % (fr * 100)
    if not nb:
        ca, cb = cal
        gb2 = tuple(statistics.fmean(gyr[i][a] for i in range(ca, cb)) for a in range(3))
        am = tuple(statistics.fmean(acc[i][a] for i in range(ca, cb)) for a in range(3))
        g_mag = M.vnorm(am)
        qs = M.solve_attitude_v2(acc, gyr, gb2, M.vunit(am), flags, dt)
        a_lin = M.to_linear(acc, qs, g_mag)

        # B：双锚点（终点静止也用上）
        db = M.zupt_two_anchor(a_lin, w0, w1, run0, run1, dt)[0]
        # C：只减一个常值残差，残差由按下前那段静止量出来（用户方案）
        la, lb = run0[0], min(run0[1], dn)
        b_hat = tuple(statistics.fmean(a_lin[i][a] for i in range(la, lb)) for a in range(3))
        a_corr = [tuple(a_lin[i][x] - b_hat[x] for x in range(3)) for i in range(n)]
        dc = M.vnorm(M.integrate_from(a_corr, w0, w1, dt)[1][-1])
        nb["bareM"] = M.vnorm(M.integrate_from(a_lin, w0, w1, dt)[1][-1])
        nb["T"] = (w1 - w0) * dt
        nb["anchorRunS"] = ((run0[1] - run0[0]) * dt, (run1[1] - run1[0]) * dt)

    return {"A": da, "B": db, "C": dc, "A_res": ra, "BC": nb,
            "holdS": (up - dn) * dt, "n": n, "dt": dt}


def main():
    names = sorted(x for x in os.listdir(CAP) if os.path.isdir(os.path.join(CAP, x)))
    print("=" * 104)
    print("实验 14：测量区间由谁决定（A 静止切段 / B 按钮窗+双锚点 / C 按钮窗+只起点）")
    print("=" * 104)
    print("%-17s %5s %6s | %9s %9s %9s | %8s %8s | %s"
          % ("采集", "真值", "按住s", "A 现役", "B 双锚点", "C 只起点", "裸积分", "C-B差", "备注"))
    print("-" * 104)

    errs = {"A": [], "B": [], "C": []}
    for nm in names:
        v = variants(os.path.join(CAP, nm))
        if v is None:
            continue
        t = TRUTH.get(nm)
        f = lambda x: "%9.2f" % (x * 100) if x is not None else "%9s" % "-"
        note = ""
        if t:
            e = [(v[k] * 100 - t) / t * 100 if v[k] is not None else float("nan")
                 for k in ("A", "B", "C")]
            for k, ee in zip(("A", "B", "C"), e):
                if v[k] is not None:
                    errs[k].append(ee)
            note = "A %+.1f%%  B %+.1f%%  C %+.1f%%" % tuple(e)
        elif v["BC"].get("why"):
            note = v["BC"]["why"]
        elif v["B"] is not None and v["C"] is not None:
            note = "裸积分 %.0f cm" % (v["BC"]["bareM"] * 100)
        mark = "★" if t else " "
        cb = ("%+.2f" % ((v["C"] - v["B"]) * 100)) if (v["B"] is not None and v["C"] is not None) else "-"
        print("%s%-16s %5s %6.2f | %s %s %s | %8s %8s | %s"
              % (mark, nm, ("%.0f" % t) if t else "-", v["holdS"],
                 f(v["A"]), f(v["B"]), f(v["C"]),
                 ("%.0f" % (v["BC"].get("bareM", 0) * 100)) if v["BC"].get("bareM") else "-",
                 cb, note))

    print("-" * 104)
    print("有真值的三条上（真值 15 cm，本身 ±3.3%）：")
    lab = {"A": "现役(静止切段)", "B": "按钮窗+双锚点", "C": "按钮窗+只起点"}
    for k in ("A", "B", "C"):
        if errs[k]:
            print("   %-16s 中位误差 %+.2f%%   最差 %+.2f%%   （%d 条）"
                  % (lab[k], statistics.median(errs[k]),
                     max(abs(e) for e in errs[k]), len(errs[k])))
    print()
    print("「C 相对 B 差多少」= 终点静止值多少钱。")


if __name__ == "__main__":
    main()
