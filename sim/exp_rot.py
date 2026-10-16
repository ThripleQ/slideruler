"""旋转对读数的影响 —— 把「感觉旋转影响很大」变成可检验的问题。

2026-10-16 夜 手操 20 cm 批次：A 组 4 条 + 夜批 22 条 = **26 条**（设备上另两条
233740 是健康检查记录、234911 被引擎拒答，非测量）。

真值与工况由用户 2026-10-17 确认：**全部 20 cm、全部平放桌面推**。所以下面的
结论在本批内是完整的（不用再挂"待确认"），但**不能外推到手持悬空推** —— 那时
起止姿态变化会大得多。

对每条算三个「旋转量」，再看误差跟谁走：
  theta_net  窗口内姿态**净**转角，由引擎返回的 q0/q1 得到（度）
  theta_int  窗口内陀螺**总**转角 ∫‖ω‖dt（度，不区分方向、含来回抖）
  leverMax   若 IMU 到测量头 r = 9 cm 且方向最不利，杠杆臂会贡献多少位移
             = 2·|r|·sin(θ_net/2)（cm）。这是**上限**，实际值取决于 r 与转轴夹角。

再试一件正事：用 26 条 + 真值 20 cm 反解 r 的三个分量。
⚠ 这是"拿真值拟合参数"，有循环论证嫌疑，只作**探索**：
  · 前提是各条 ΔR 的方向足够分散 —— MEMORY 记着「ΔR 对任何单一旋转必然奇异」，
    ΔR·r = (Rot(n,θ) − I)·R₀r，当 R₀r ∥ n 时恒为 0，那个方向就是零空间。
  · 所以先看 AᵀA 的谱，条件数不够就直接说"解不出来"，不许硬给个数。
  · 真值自身 ±1 cm，判据是"解出来的 r 是否落在物理合理范围"，不是小数点后几位。
"""
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure2
from exp_real import load, ROOT, HAND20
from pipeline import qconj, qmul, qrot, vnorm, vunit

TRUTH_CM = 20.0                      # 用户手操、自报 ±1 cm
TRUTH = TRUTH_CM / 100.0
LEVER_PROBE_M = 0.09                 # 只用来估算杠杆臂影响的量级

# 真值批次登记在 exp_real.HAND20 —— 这里不另立列表，免得多头登记漂移。
# 用户 2026-10-17 确认：全部 20 cm、**全部平放桌面推**。
CASES = sorted(HAND20)


def rot_cols(q):
    """R(q) 的三列：R·e_j = qrot(q, e_j)。qrot 是 body→world，所以列就是世界系基像。"""
    return [qrot(q, e) for e in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))]


def d_rot_mat(q0, q1):
    """ΔR = R(q1) − R(q0)，按列存：dR[j] 是第 j 列。"""
    c0, c1 = rot_cols(q0), rot_cols(q1)
    return [tuple(c1[j][i] - c0[j][i] for i in range(3)) for j in range(3)]


def matvec_by_cols(cols, v):
    return tuple(sum(cols[j][i] * v[j] for j in range(3)) for i in range(3))


def net_angle_deg(q0, q1):
    """q0→q1 的转角。dq = q1·q0⁻¹，θ = 2·acos|w|（与坐标系约定无关）。"""
    dq = qmul(q1, qconj(q0))
    return math.degrees(2.0 * math.acos(min(1.0, abs(dq[0]))))


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx < 1e-12 or sy < 1e-12:
        return float("nan")
    return sum((xs[i] - mx) * (ys[i] - my) for i in range(n)) / (sx * sy)


def solve3(A, b):
    """最小二乘 x=(AᵀA)⁻¹Aᵀb。A 是行列表（每行长度 3）。返回 (x, ATA)。"""
    n = len(A)
    ATA = [[sum(A[k][i] * A[k][j] for k in range(n)) for j in range(3)] for i in range(3)]
    ATb = [sum(A[k][i] * b[k] for k in range(n)) for i in range(3)]
    M = [ATA[i][:] + [ATb[i]] for i in range(3)]
    for c in range(3):                                   # 带部分主元的高斯消元
        p = max(range(c, 3), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-15:
            return None, ATA
        M[c], M[p] = M[p], M[c]
        for r in range(3):
            if r != c:
                f = M[r][c] / M[c][c]
                for k in range(c, 4):
                    M[r][k] -= f * M[c][k]
    return [M[i][3] / M[i][i] for i in range(3)], ATA


def eig_sym3(M):
    """3×3 对称阵特征值（雅可比迭代，够用）。"""
    a = [row[:] for row in M]
    for _ in range(60):
        off = abs(a[0][1]) + abs(a[0][2]) + abs(a[1][2])
        if off < 1e-15:
            break
        for p, q in ((0, 1), (0, 2), (1, 2)):
            if abs(a[p][q]) < 1e-18:
                continue
            theta = 0.5 * math.atan2(2 * a[p][q], a[q][q] - a[p][p])
            c, s = math.cos(theta), math.sin(theta)
            for k in range(3):
                akp, akq = a[k][p], a[k][q]
                a[k][p] = c * akp - s * akq
                a[k][q] = s * akp + c * akq
            for k in range(3):
                apk, aqk = a[p][k], a[q][k]
                a[p][k] = c * apk - s * aqk
                a[q][k] = s * apk + c * aqk
    return sorted((a[i][i] for i in range(3)), reverse=True)


def main():
    rows = []
    for nm in CASES:
        got = load(os.path.join(ROOT, nm))
        if not got:
            print("跳过（load 失败）:", nm)
            continue
        acc, gyr, dn, up, dt = got
        m = measure2.measure(acc, gyr, dn, up, dt=dt, merge_mode="press",
                             head_r=None, gate_motion=True)
        if not m.get("ok"):
            print("拒答:", nm, m.get("reason"))
            continue
        w0, w1 = m["window"]
        q0, q1 = m["q0"], m["q1"]
        th_net = net_angle_deg(q0, q1)
        th_int = math.degrees(sum(vnorm(gyr[i]) for i in range(w0, w1)) * dt)
        lever_max = 2.0 * LEVER_PROBE_M * math.sin(math.radians(th_net) / 2.0) * 100.0
        p_imu = m["pImuVec"]
        rows.append({
            "nm": nm, "theta_net": th_net, "theta_int": th_int,
            "lever_max": lever_max, "read": m["distanceM"] * 100,
            "err": m["distanceM"] * 100 - TRUTH_CM,
            "bare": m["bareM"] * 100, "gap": m["modelGapM"] * 100,
            "drift": m["anchorDrift"], "motionS": m["motionS"],
            "nseg": m["nMerged"], "p_imu": p_imu, "norm_p": vnorm(p_imu) * 100,
            "dR": d_rot_mat(q0, q1), "plateau": m["plateau"],
        })

    print("旋转 vs 误差（手操 20 cm，%d 条）" % len(rows))
    print("=" * 100)
    print("%-16s %7s %7s %8s %8s %7s %7s %7s %6s"
          % ("采集", "净转角°", "总转角°", "杠杆上限", "读数cm", "误差cm", "裸积cm", "锚漂", "段"))
    print("-" * 100)
    for r in sorted(rows, key=lambda x: x["theta_net"]):
        print("%-16s %7.1f %7.1f %8.2f %8.2f %+7.2f %7.2f %7.4f %6d"
              % (r["nm"], r["theta_net"], r["theta_int"], r["lever_max"],
                 r["read"], r["err"], r["bare"], r["drift"], r["nseg"]))

    print("-" * 100)
    abs_err = [abs(r["err"]) for r in rows]
    reads = [r["read"] for r in rows]
    print("读数 均值 %.2f 中位 %.2f  σ %.2f  最小 %.2f  最大 %.2f cm"
          % (statistics.fmean(reads), statistics.median(reads), statistics.pstdev(reads),
             min(reads), max(reads)))
    print("误差 均值 %+.2f  中位 %+.2f  |误差| 均值 %.2f  RMS %.2f cm"
          % (statistics.fmean([r["err"] for r in rows]),
             statistics.median([r["err"] for r in rows]),
             statistics.fmean(abs_err),
             math.sqrt(statistics.fmean([e * e for e in abs_err]))))

    th = [r["theta_net"] for r in rows]
    ti = [r["theta_int"] for r in rows]
    lv = [r["lever_max"] for r in rows]
    gp = [r["gap"] for r in rows]
    dr = [abs(r["drift"]) for r in rows]
    print("\n相关（Pearson）：")
    print("  |误差| ~ 净转角      %.3f" % pearson(th, abs_err))
    print("  |误差| ~ 总转角      %.3f" % pearson(ti, abs_err))
    print("  |误差| ~ 杠杆臂上限  %.3f" % pearson(lv, abs_err))
    print("  |误差| ~ |锚点漂移|  %.3f" % pearson(dr, abs_err))
    print("  |误差| ~ 修正量gap   %.3f" % pearson(gp, abs_err))
    print("  |误差| ~ 运动时长    %.3f" % pearson([r["motionS"] for r in rows], abs_err))
    print("  净转角 ~ 总转角      %.3f" % pearson(th, ti))

    if rows:
        half = max(1, len(rows) // 2)
        srt = sorted(rows, key=lambda x: x["theta_net"])
        lo, hi = srt[:half], srt[len(srt) - half:]
        rms = lambda sel: math.sqrt(statistics.fmean([x["err"] ** 2 for x in sel]))
        print("\n按净转角二分：低半 (%s~%.1f°) |误差| RMS %.2f cm  ｜ 高半 (%.1f~%s°) |误差| RMS %.2f cm"
              % ("%.1f" % lo[0]["theta_net"], lo[-1]["theta_net"], rms(lo),
                 hi[0]["theta_net"], "%.1f" % hi[-1]["theta_net"], rms(hi)))

    print("\n" + "=" * 100)
    print("反解杠杆臂 r（用真值 20 cm）")
    print("-" * 100)
    A, b, used = [], [], []
    for r in rows:
        u = vunit(r["p_imu"])
        dR = r["dR"]
        A.append([sum(u[k] * dR[j][k] for k in range(3)) for j in range(3)])
        b.append(TRUTH - r["norm_p"] / 100.0)
        used.append(r["nm"])
    x, ATA = solve3(A, b)
    if x is None:
        print("AᵀA 奇异，解不出来 —— 这批数据的 ΔR 方向不足以定 r。")
    else:
        ev = eig_sym3(ATA)
        cond = ev[0] / ev[2] if ev[2] > 1e-18 else float("inf")
        print("ΔR 方向多样性：AᵀA 特征值 %.4g / %.4g / %.4g   条件数 %.1f"
              % (ev[0], ev[1], ev[2], cond))
        print("解出 r = (%.2f, %.2f, %.2f) cm   模长 %.2f cm"
              % (x[0] * 100, x[1] * 100, x[2] * 100, vnorm(x) * 100))
        resid = [1.0]  # placeholder
        print("（r 是 body 系里 IMU→测量头的向量；物理量级预期几 cm，超过 20 cm 就要怀疑）")
        worst = max(range(len(A)), key=lambda k: abs(b[k]))
        print("参与拟合最大偏差 %.2f cm（%s）—— 真值本身 ±1 cm" % (abs(b[worst]) * 100, used[worst]))
        # 线性化合理性自检
        print("线性化自检：|ΔR·r| 最大 %.2f cm vs ‖p_imu‖ 中位 %.2f cm"
              % (max(vnorm(matvec_by_cols(r["dR"], x)) for r in rows) * 100,
                 statistics.median([r["norm_p"] for r in rows])))
    print("\n⚠ 这一步是拿真值拟合参数，**只作探索**，不能当结论 —— 见文件头说明。")


if __name__ == "__main__":
    main()
