"""实验 20：终点锚点到底停没停？

## 起因

20261011_143734（真值 15 cm）原来报 1.44 cm。诊断（sim/diag_anchor.py）显示
**裸积分 16.19 cm 本来是对的**，是双锚点扣掉了 15.42 cm —— 而扣掉的那笔不是漂移：

    终点锚点段（0.143 s）内积分速度稳定在 −0.31 m/s，运动段速度峰值才 0.314 m/s。

也就是说，标注为「静止」的那一段里，手机正以接近滑动的速度**匀速**移动。
「松手」不等于「停住」，而**匀速移动与静止在加速度计上不可分辨**（测的是比力，
匀速时比力就是重力），静止检测必然把它判成静止。v3 那条教训（匀速段把一次滑动
切碎）当时只用来改**运动段**判据，没意识到它同样毁掉**终点锚点** —— 而终点锚点
正是双锚点里"把漂移钉回零"的那一半。同一条的起点锚点 |V0| = 0.007 m/s（真静止），
问题只在终点一侧。

## 试过的判据：用 att_budget 当分母 —— **被否**

位移上界 D = att_budget，对应恒定残差加速度 A_sup = 2D/T²，速度上界 V_sup = 2D/T。
「锚点里要去掉的那点速度不可能超过 V_sup」听起来有物理依据，但实测不成立：

    143734  比值 2.84  触发（对的）
    100529  比值 1.54  触发（**误伤** —— 它是 +1.6% 的正确读数）
    122239  比值 0.08  漏掉（锚点速度 1.02 m/s，明显在动）
    122244  比值 0.50  漏掉（锚点速度 13.5 m/s！）

原因：V_sup ∝ σ_bg，而 σ_bg 由静止段切 4 段、仅 3 个差值估计，本身方差极大 ——
全库 V_sup 跨 4 个数量级，等于没有判别力。**不能把它当拒答门槛。**

## 采用的诊断：锚点速度 / 运动段速度峰值

    anchorSuspect := v_pk > 0.05 m/s 且 |V_anchor| > 0.5·v_pk

比值 ≈ 1 的含义很直白：标为"静止"的那段里，速度等于整个运动过程的**最大**速度 ——
手机在窗口末端仍在以峰值速度移动，完全没减速。实测分界干净：

    有真值三条正常   100444 0.09 / 100459 0.16 / 100529 0.06
    有真值三条可疑   143719 0.98 / 143725 0.90 / 143734 1.01

**但它只诊断、不拒答。** 143725 的 0.90 同样可疑，却给出 +5% 的正确读数
（两个错误互相抵消）；拿 3~6 条真值去拟合一条分界线，等于用真值调参，本项目不做。

真正被拒答的是 **modelSensitive**（见 measure2.py §6e）—— 那个判据本来就在说
"这个数由模型决定、不是由数据决定"，那它就不该给数。143734 落进那里。

## 追加（2026-10-16）：att_budget 本身也退役了

上面那条"被否的判据"的分母 `att_budget` 已从引擎删除，理由与它作为判据无关：
① 前提被证伪 —— 库里 `20261011_085650` 是 76 s 纯静止录制，那上面零偏不游走
   （5 s 窗口均值全在 ±3e-5 rad/s 内、Allan 从 10 ms 到 9.7 s 单调下降、无平台）；
② 它算出的预算与实际误差**反相关**（预算 5.4 cm 的 143734 错得最狠 −90%，
   预算 35.7 cm 的 143725 反而 +5%），一个被自己的数据违反 2.5 倍的"上界"不是上界。
引擎里现在只剩 `startTremorRadS / startTremorTauS`（起手前静止段的角抖动）。

本脚本为了**留痕**，仍按已退役的公式现场重算 V_sup，好让下表逐字可复现：
    D = (1/6)·g·σ·√(T/τ)·T³   ⟹   V_sup = 2D/T = (1/3)·g·σ·√(T/τ)·T²
表里那一列标注为"已退役"。**不要**再拿它做任何判据。

用法：python sim/exp_anchor2.py
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2 as M
from exp_real import load, TRUTH_CM, tier_of, ROOT


def main():
    names = sorted(x for x in os.listdir(ROOT) if os.path.isdir(os.path.join(ROOT, x)))
    print("%-16s %8s %8s %8s %8s %8s  %s" %
          ("采集", "锚点速度", "运动峰值", "比值", "V_sup退役", "被否比值", "判定 / 真值"))
    print("-" * 112)
    n_susp = n_tot = n_rej_old = 0
    for nm in names:
        d = load(os.path.join(ROOT, nm))
        if d is None:
            continue
        acc, gyr, dn, up, dt = d
        m = M.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                      head_r=(0.0, -0.09, 0.0))
        ad = m.get("anchorDriftMps", m.get("anchorDrift"))
        vp = m.get("vPeakMps")
        if ad is None or vp is None:
            continue
        n_tot += 1
        ratio = ad / vp if vp > 1e-9 else 0.0
        T = m.get("T") or m.get("windowS") or 0.0
        # 已退役的 att_budget 公式（见文件头"追加"），只为留痕、逐字复现旧表：
        #   D = (1/6)·g·σ·√(T/τ)·T³  ⟹  V_sup = 2D/T
        sig = m.get("startTremorRadS") or 0.0
        tau = m.get("startTremorTauS") or 0.0
        g = m.get("gMag") or 0.0
        D = (g * sig * math.sqrt(T / tau) * T ** 3 / 6.0) if (sig and tau and T > 0) else 0.0
        sup = 2.0 * D / T if T > 0 else 0.0
        rej_ratio = ad / sup if sup > 1e-12 else float("inf")
        susp = m.get("anchorSuspect", False)
        n_susp += 1 if susp else 0
        n_rej_old += 1 if rej_ratio > 1.0 else 0
        t = TRUTH_CM.get(nm)
        note = ("可疑" if susp else "")
        if t:
            note += ("  " if note else "") + "真值%.1f(%s)" % (t, tier_of(nm))
            if m.get("ok"):
                note += " %+.1f%%" % ((m["distanceM"] * 100 - t) / t * 100)
            else:
                note += " 已拒答"
        print("%-16s %8.4f %8.4f %8.2f %8.4f %8.2f  %s"
              % (nm, ad, vp, ratio, sup, rej_ratio, note))
    print("-" * 112)
    print("可评估 %d 条：锚点可疑 %d 条（只诊断，不拒答）；"
          "若改用 att_budget 版判据（倒数第二列 >1）会拒 %d 条 —— 已否，见文件头"
          % (n_tot, n_susp, n_rej_old))


if __name__ == "__main__":
    main()
