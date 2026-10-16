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

用法：python sim/exp_anchor2.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2 as M
from exp_real import load, TRUTH_CM, HARD, ROOT


def main():
    names = sorted(x for x in os.listdir(ROOT) if os.path.isdir(os.path.join(ROOT, x)))
    print("%-16s %8s %8s %8s %8s %8s  %s" %
          ("采集", "锚点速度", "运动峰值", "比值", "V_sup", "被否比值", "判定 / 真值"))
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
        sup = 2.0 * m.get("attBudgetM", 0.0) / T if T > 0 else 0.0
        rej_ratio = ad / sup if sup > 1e-12 else float("inf")
        susp = m.get("anchorSuspect", False)
        n_susp += 1 if susp else 0
        n_rej_old += 1 if rej_ratio > 1.0 else 0
        t = TRUTH_CM.get(nm)
        note = ("可疑" if susp else "")
        if t:
            note += ("  " if note else "") + "真值%.1f%s" % (t, "(硬)" if nm in HARD else "(估)")
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
