package com.thripleq.slideruler.probe

import kotlin.math.abs
import kotlin.math.sqrt

/** 单条流的统计量。 */
class StreamStats(
    val n: Int,
    val mean: List<Double>,
    val std: List<Double>,
    val normMean: Double,
    val normStd: Double,
    val rateHz: Double,
    val spanS: Double,
    val dtMedMs: Double,
    val dtP05Ms: Double,
    val dtP95Ms: Double,
    val dtMaxMs: Double,
    val drops: Int,
    val hist: List<Int>,
)

/**
 * 体检 A（传感器静态体检）—— 对应 docs/implementation-guide.md §2「怎么验证」四条：
 *   ① ‖accel‖ 应 ≈ 9.81（偏差 > 0.3 说明这颗传感器有问题）
 *   ② 三轴均值 = 该机的零偏水平（记录下来，机型档案）
 *   ③ 三轴标准差 = 噪声水平（0.01~0.05 m/s² 都正常，> 0.08 这台机器难做）
 *   ④ timestamp 间隔直方图：应集中在 10 ms，尾巴长说明调度不稳
 * 外加两条我们自己关心的：实测采样率、丢样次数（dt > 50 ms）。
 */
object HealthCheck {

    private const val TARGET_HZ = 100.0
    private const val G = 9.81

    fun analyse(samples: List<Sample>): StreamStats {
        val n = samples.size
        require(n >= 2) { "样本不足" }

        val mean = DoubleArray(3)
        for (s in samples) {
            mean[0] += s.x
            mean[1] += s.y
            mean[2] += s.z
        }
        for (a in 0..2) mean[a] /= n

        val varSum = DoubleArray(3)
        var normSum = 0.0
        for (s in samples) {
            val dx = s.x - mean[0]
            val dy = s.y - mean[1]
            val dz = s.z - mean[2]
            varSum[0] += dx * dx
            varSum[1] += dy * dy
            varSum[2] += dz * dz
            normSum += norm(s)
        }
        val std = List(3) { sqrt(varSum[it] / (n - 1)) }
        val normMean = normSum / n

        var nv = 0.0
        for (s in samples) {
            nv += (norm(s) - normMean) * (norm(s) - normMean)
        }
        val normStd = sqrt(nv / (n - 1))

        val dts = DoubleArray(n - 1)
        for (i in 1 until n) dts[i - 1] = (samples[i].tNs - samples[i - 1].tNs) / 1e6
        val sorted = dts.sortedArray()

        val spanS = (samples[n - 1].tNs - samples[0].tNs) / 1e9
        val rate = if (spanS > 0) (n - 1) / spanS else 0.0

        val hist = IntArray(5)
        for (d in dts) hist[bucket(d)]++

        return StreamStats(
            n = n,
            mean = mean.toList(),
            std = std,
            normMean = normMean,
            normStd = normStd,
            rateHz = rate,
            spanS = spanS,
            dtMedMs = percentile(sorted, 0.50),
            dtP05Ms = percentile(sorted, 0.05),
            dtP95Ms = percentile(sorted, 0.95),
            dtMaxMs = sorted.last(),
            drops = dts.count { it > ImuRecorder.DROP_MS },
            hist = hist.toList(),
        )
    }

    private fun norm(s: Sample): Double {
        val x = s.x.toDouble()
        val y = s.y.toDouble()
        val z = s.z.toDouble()
        return sqrt(x * x + y * y + z * z)
    }

    private fun bucket(dtMs: Double): Int = when {
        dtMs < 8.0 -> 0
        dtMs < 12.0 -> 1
        dtMs < 20.0 -> 2
        dtMs < 50.0 -> 3
        else -> 4
    }

    private fun percentile(sorted: DoubleArray, q: Double): Double {
        if (sorted.isEmpty()) return 0.0
        val idx = ((sorted.size - 1) * q).toInt().coerceIn(0, sorted.size - 1)
        return sorted[idx]
    }

    /** 生成给人看的体检报告（纯文本，可直接贴进机型档案）。 */
    fun buildReport(acc: StreamStats, gyr: StreamStats, mode: String): String = buildString {
        appendLine("体检 · $mode")
        appendLine("acc ${acc.n} 样本 / gyr ${gyr.n} 样本 · 时长 ${fmt(acc.spanS, 2)} s")
        appendLine()
        appendLine("── 加速度计（m/s²）──")
        val dev = abs(acc.normMean - G)
        appendLine("  ‖a‖ 均值  ${fmt(acc.normMean, 4)}   偏离 9.81  ${fmt(dev, 4)}   ${mark(dev < 0.10, dev < 0.30)}")
        appendLine("  零偏   x ${sgn(acc.mean[0])}  y ${sgn(acc.mean[1])}  z ${sgn(acc.mean[2])}")
        appendLine("  噪声σ  x ${fmt(acc.std[0], 4)}  y ${fmt(acc.std[1], 4)}  z ${fmt(acc.std[2], 4)}   ${noiseMark(acc.std.max())}")
        appendLine("  实测采样率 ${fmt(acc.rateHz, 1)} Hz   间隔中位 ${fmt(acc.dtMedMs, 2)} ms"
            + "  P05 ${fmt(acc.dtP05Ms, 2)}  P95 ${fmt(acc.dtP95Ms, 2)}  max ${fmt(acc.dtMaxMs, 1)}")
        appendLine("  丢样(>${ImuRecorder.DROP_MS.toInt()}ms) ${acc.drops} 次   ${histogram(acc)}")
        appendLine("  ${rateMark(acc)}")
        appendLine()
        appendLine("── 陀螺仪（rad/s）──")
        appendLine("  零偏   x ${sgn(gyr.mean[0])}  y ${sgn(gyr.mean[1])}  z ${sgn(gyr.mean[2])}")
        appendLine("         = ${fmt(Math.toDegrees(gyr.mean[0]), 3)} / ${fmt(Math.toDegrees(gyr.mean[1]), 3)}"
            + " / ${fmt(Math.toDegrees(gyr.mean[2]), 3)} °/s")
        appendLine("  噪声σ  x ${fmt(gyr.std[0], 5)}  y ${fmt(gyr.std[1], 5)}  z ${fmt(gyr.std[2], 5)}")
        appendLine("  范数均值 ${fmt(gyr.normMean, 5)} rad/s（静止时应 ≈ 0，非 0 部分就是零偏水平）")
        appendLine("  实测采样率 ${fmt(gyr.rateHz, 1)} Hz   间隔中位 ${fmt(gyr.dtMedMs, 2)} ms"
            + "  P95 ${fmt(gyr.dtP95Ms, 2)}  max ${fmt(gyr.dtMaxMs, 1)}   丢样 ${gyr.drops} 次")
        appendLine()
        appendLine("── 判读 ──")
        appendLine("  ${verdict(acc, gyr)}")
    }

    private fun histogram(s: StreamStats): String =
        "直方图 <8:${s.hist[0]}  8–12:${s.hist[1]}  12–20:${s.hist[2]}  20–50:${s.hist[3]}  >50:${s.hist[4]}"

    private fun rateMark(s: StreamStats): String = when {
        s.drops > 0 -> "✗ 出现 ${s.drops} 次丢样：算法层必须按实测 dt 积分并标记这些段"
        s.dtP95Ms > 13.0 -> "△ P95 间隔 ${fmt(s.dtP95Ms, 1)} ms，调度不稳；仍必须按实测 dt 积分"
        else -> "✓ 采样节奏稳定"
    }

    private fun noiseMark(maxStd: Double): String = when {
        maxStd <= 0.05 -> "✓ 正常（指南说 0.01~0.05 都正常）"
        maxStd <= 0.08 -> "△ 偏大，这台机器会更难"
        else -> "✗ 过大（>0.08），这台机器难做"
    }

    private fun mark(ok: Boolean, warn: Boolean): String = when {
        ok -> "✓"
        warn -> "△"
        else -> "✗ 这颗传感器有问题（偏差 > 0.3）"
    }

    private fun verdict(acc: StreamStats, gyr: StreamStats): String {
        val dev = abs(acc.normMean - G)
        val noise = acc.std.max()
        val rateOk = abs(acc.rateHz - TARGET_HZ) < 5.0
        val bad = dev > 0.30 || noise > 0.08 || !rateOk || acc.drops > 0
        return when {
            !bad -> "本机传感器合格：文档预算（噪声 0.02 m/s²）适用于这台机器，可以进入下一阶段。"
            else -> "本机存在异常项（见上），先别采滑动数据——要先拿到「这台机器的真实数字」再继续。"
        }
    }

    private fun fmt(v: Double, digits: Int): String {
        // 不用 "%.Nf".format：它受默认 Locale 影响，某些区域会把小数点写成逗号
        val s = "%.${digits}f".format(java.util.Locale.US, v)
        return s
    }

    private fun sgn(v: Double): String {
        val s = fmt(v, 4)
        return if (v >= 0) "+$s" else s
    }
}
