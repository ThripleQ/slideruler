"""实验 18：真机静止段的噪声结构 vs 仿真器 —— 给仿真器的误差模型做标定。

## 为什么必须做

`exp_zupt.py` 里同一个改动（残差模型 线性→二次）在两侧给出**相反**判决：
仿真说"二次好很多"（T=1.4 s：+23%→+0%），真机说"二次差很多"（143719：−21%→−78%）。
说明 `SimCfg.gyr_bias_rw` 那个**白噪声随机游走**不是 K40 真实零偏误差的形状。
在把仿真器的误差模型标定到与真机一致之前，仿真**不能**用来 accept/reject 估计器。

## 怎么做

Allan deviation（重叠式，短序列也稳）：把静止段按 τ 个样本分组，组内均值之差的
RMS/√2 就是该 τ 上的 σ。三种典型形状的斜率完全不同：

    σ(τ) ∝ τ^(-1/2)   白噪声（角度随机游走 / ARW）      —— 每样本独立噪声
    σ(τ) ∝ τ^0        **偏置不稳定**（1/f，有一个平台）  —— 抹不掉、平均不掉
    σ(τ) ∝ τ^(+1/2)   速率随机游走（零偏在游走）         —— 分段均值会漂

本项目的 T³ 定律里，主宰项就是"偏置不稳定 → 零偏游走"，所以要看的是
**中长 τ 上有没有平台、平台多高**，而不是短 τ 的白噪声大小。

输出两张表：陀螺、加速度。真机与仿真并排，同一个 τ 网格。
"""
import math
import os
import statistics
import sys

sys.path.insert(0, '.')

import measure2 as M
from exp_real import ROOT


def static_segments(series, min_len=40):
    """从一条采集的静止标记里取出所有静止段的样本下标。"""
    return [(a, b) for (a, b) in series if b - a >= min_len]


def allan_curve(rows, taus=(1, 2, 4, 8, 16, 32, 64), min_pairs=3, robust=False):
    """重叠式 Allan deviation。

    rows：[[(x,y,z), ...], ...] 若干条独立的静止段（每条段内不跨接）。
    taus：τ 的网格（单位：样本数）。
    返回 {τ(样本数): σ}，σ 按三轴合并。

    robust=False（默认）：σ = sqrt(mean(d²)/2) —— 经典 Allan，对离群敏感。
    robust=True：σ = median(|d|/√2) —— 抗离群。本项目**必须**看这一条：
      真机静止段里有假静止（匀速段被静止检测误判），它们把 RMS 撑高一个量级。
      两者之差 = "上升段到底是真的零偏游走，还是少数几条在动的段撑起来的"。
    """
    out = {}
    for tau in taus:
        sq, ab = [], []
        for seg in rows:
            n = len(seg)
            if n < tau * min_pairs:
                continue
            for x in range(3):
                cum = [0.0]
                for v in seg:
                    cum.append(cum[-1] + v[x])
                for i in range(0, n - 2 * tau + 1):
                    a = (cum[i + tau] - cum[i]) / tau
                    b = (cum[i + 2 * tau] - cum[i + tau]) / tau
                    sq.append((a - b) ** 2)
                    ab.append(abs(a - b) / math.sqrt(2.0))
        if sq:
            out[tau] = (statistics.median(ab) if robust
                        else math.sqrt(statistics.fmean(sq) / 2.0))
    return out


def collect_real():
    """所有真机采集里的静止段，分陀螺/加速度两组。"""
    gyr_rows, acc_rows = [], []
    n_files = 0
    total_s = 0.0
    for nm in sorted(os.listdir(ROOT)):
        d = os.path.join(ROOT, nm)
        if not os.path.isdir(d):
            continue
        try:
            acc_ts, acc = M.read_pairs(os.path.join(d, "acc.csv"))
            gyr_ts, gyr = M.read_pairs(os.path.join(d, "gyr.csv"))
        except Exception:
            continue
        if len(acc) < 200:
            continue
        n_files += 1
        # 静止检测用最小依赖的版本：acc 局部标准差 + 陀螺模长（与 measure2 同源）
        n = min(len(acc), len(gyr))
        acc, gyr = acc[:n], gyr[:n]
        ts = acc_ts[:n]
        dt = (ts[-1] - ts[0]) / 1e9 / (len(ts) - 1)
        flags, *_ = M.static_flags(acc, gyr, (0.0, 0.0, 0.0), dt)
        segs = static_segments([(a, b) for a, b in _runs(flags)], min_len=40)
        for (a, b) in segs:
            # 去掉段均值（Allan 只看变化，均值差会污染）
            for src, dst in ((gyr, gyr_rows), (acc, acc_rows)):
                s = src[a:b]
                mu = [statistics.fmean(p[x] for p in s) for x in range(3)]
                dst.append([(p[0] - mu[0], p[1] - mu[1], p[2] - mu[2]) for p in s])
            total_s += (b - a) * dt
    return gyr_rows, acc_rows, n_files, total_s


def _runs(flags):
    a = None
    for i, f in enumerate(flags):
        if f and a is None:
            a = i
        elif not f and a is not None:
            yield (a, i)
            a = None
    if a is not None:
        yield (a, len(flags))


def collect_sim(cfg_kw, n_seg=40, seg_len=140, hz=105.1, seed=1):
    """同长度量级的仿真静止段（段长 140 样本 ≈ 真机最长静止段）。"""
    import random
    from simgen import SimCfg
    cfg = SimCfg(hz=hz, **cfg_kw)
    rng = random.Random(seed)
    gyr_rows, acc_rows = [], []
    ba = [0.0, 0.0, 0.0]            # 随机游走分量
    bg = [0.0, 0.0, 0.0]            # 有界（Gauss-Markov）分量
    sg, tc = getattr(cfg, "gyr_bias_gm", (0.0, 0.0))
    dt = 1.0 / cfg.hz
    a_gm = math.exp(-dt / tc) if tc > 0 else 0.0
    q_gm = sg * math.sqrt(max(0.0, 1.0 - a_gm * a_gm)) if tc > 0 else 0.0
    for k in range(n_seg):
        g, a = [], []
        for i in range(seg_len):
            for x in range(3):
                ba[x] += cfg.gyr_bias_rw * math.sqrt(dt) * rng.gauss(0, 1)
                bg[x] = a_gm * bg[x] + q_gm * rng.gauss(0, 1)
            g.append(tuple(ba[x] + bg[x] + cfg.gyr_bias[x]
                           + cfg.gyr_noise * rng.gauss(0, 1) for x in range(3)))
            a.append(tuple(cfg.acc_bias[x] + cfg.acc_noise * rng.gauss(0, 1)
                           for x in range(3)))
        for src, dst in ((g, gyr_rows), (a, acc_rows)):
            mu = [statistics.fmean(p[x] for p in src) for x in range(3)]
            dst.append([(p[0] - mu[0], p[1] - mu[1], p[2] - mu[2]) for p in src])
    return gyr_rows, acc_rows


def show(title, curve, unit, dt, curve_med=None):
    """打印 Allan 曲线。给了 curve_med 就把 RMS 与中位数版本并排，并给比值。

    比值列是**判别力所在**：RMS/中位数 ≈1 说明分布对称、没有离群段；
    远大于 1（真机陀螺在 τ≤76 ms 上是 1.4~2.6）说明这一段的"上升"
    主要是少数在动的假静止段撑起来的，不是零偏本身在游走。
    """
    print("  %s （σ 单位 %s）" % (title, unit))
    ts = list(curve)
    print("    τ(ms)  " + "".join("%10s" % ("%.0f" % (t * dt * 1000)) for t in ts))
    print("    σ(RMS) " + "".join("%10.6g" % curve[t] for t in ts))
    print("    斜率   " + "".join(
        "%10s" % ("—" if i == 0 else "%.2f" % (
            (math.log(curve[ts[i]]) - math.log(curve[ts[i - 1]]))
            / (math.log(ts[i]) - math.log(ts[i - 1]))))
        for i in range(len(ts))))
    if curve_med is not None:
        print("    σ(中位)" + "".join("%10.6g" % curve_med[t] for t in ts))
        print("    斜率   " + "".join(
            "%10s" % ("—" if i == 0 else "%.2f" % (
                (math.log(curve_med[ts[i]]) - math.log(curve_med[ts[i - 1]]))
                / (math.log(ts[i]) - math.log(ts[i - 1]))))
            for i in range(len(ts))))
        print("    RMS/中位" + "".join("%10.2f" % (curve[t] / curve_med[t])
                                     for t in ts))
    if len(ts) >= 2:
        pk = max(ts, key=lambda t: curve[t])
        print("    峰值 %.5g @ τ=%.0f ms  （白噪声 −0.5 / 有界游走：先升后降 / 无界游走 +0.5）"
              % (curve[pk], pk * dt * 1000))
        if curve_med is not None:
            pkm = max(ts, key=lambda t: curve_med[t])
            print("    峰值 %.5g @ τ=%.0f ms（中位数版）"
                  % (curve_med[pkm], pkm * dt * 1000))


def longest_static(ROOT_dir):
    """全库最长的一条连续静止段。返回 (名字, a, b, acc, gyr, dt, 总样本)。"""
    best = (0, None, None)
    for nm in sorted(os.listdir(ROOT_dir)):
        d = os.path.join(ROOT_dir, nm)
        if not os.path.isdir(d):
            continue
        try:
            acc, gyr, dn, up, dt = __import__("exp_real").load(d)
        except Exception:
            continue
        n = min(len(acc), len(gyr))
        acc, gyr = acc[:n], gyr[:n]
        flags, *_ = M.static_flags(acc, gyr, (0.0, 0.0, 0.0), dt)
        a = None
        for i, f in enumerate(list(flags) + [False]):
            if f and a is None:
                a = i
            elif not f and a is not None:
                if i - a > best[0]:
                    best = (i - a, nm, (a, i))
                a = None
    ln, nm, (a, b) = best[0], best[1], best[2]
    acc, gyr, dn, up, dt = __import__("exp_real").load(os.path.join(ROOT_dir, nm))
    return nm, a, b, acc, gyr, dt, len(acc)


def main():
    print("=" * 88)
    print("真机（全部采集的静止段）")
    gyr_rows, acc_rows, nf, total_s = collect_real()
    print("  采集 %d 条，静止段 %d 段 / 加速度 %d 段，累计静止 %.1f s"
          % (nf, len(gyr_rows), len(acc_rows), total_s))
    dt = 0.009514
    show("陀螺", allan_curve(gyr_rows), "rad/s", dt,
         allan_curve(gyr_rows, robust=True))
    show("加速度", allan_curve(acc_rows), "m/s²", dt,
         allan_curve(acc_rows, robust=True))

    print()
    print("=" * 88)
    nm, a, b, acc, gyr, dtx, n_all = longest_static(ROOT)
    print("最长单条静止段：%s  样本 [%d,%d) = %d 样本 = %.1f s（整条录制 %.1f s）"
          % (nm, a, b, b - a, (b - a) * dtx, n_all * dtx))
    print("  —— 这一段直接回答\"零偏 1~60 s 尺度有界还是无界\"，不必等新的长静止采集")
    seg = [gyr[i] for i in range(a, b)]
    taus = (1, 4, 16, 64, 256, 1024, 2048)
    show("陀螺（整段，未扣均值）", allan_curve([seg], taus=taus), "rad/s", dtx,
         allan_curve([seg], taus=taus, robust=True))

    print()
    print("=" * 88)
    print("仿真器（gyr_noise=0.002）")
    for tag, kw in (
        ("无游走（只有白噪声）", dict(gyr_bias_rw=0.0)),
        ("无界随机游走 rw=0.02（当前在用）", dict(gyr_bias_rw=0.02)),
        ("有界游走 GM(0.0045 rad/s, 50 ms)", dict(gyr_bias_gm=(0.0045, 0.05))),
        ("有界 + 小游走 rw=0.003", dict(gyr_bias_rw=0.003, gyr_bias_gm=(0.0045, 0.05))),
        ("有界 三倍幅值 GM(0.0135, 50 ms) + rw=0.003", dict(
            gyr_bias_rw=0.003, gyr_bias_gm=(0.0135, 0.05))),
    ):
        gyr_rows_s, acc_rows_s = collect_sim(dict(
            gyr_noise=0.002, acc_noise=0.008, gyr_bias=(0, 0, 0), **kw))
        print("  --- %s ---" % tag)
        show("陀螺", allan_curve(gyr_rows_s), "rad/s", 1 / 105.1,
             allan_curve(gyr_rows_s, robust=True))


if __name__ == "__main__":
    main()
