package com.thripleq.slideruler.probe

import android.content.Context
import android.os.Build
import java.io.BufferedWriter
import java.io.File
import java.text.SimpleDateFormat
import java.util.Date
import java.util.Locale

/**
 * 一次触摸标记（按住式测量的边界）。
 *
 * 两个时间基准必须同时留档，因为它们是**不同**的时钟：
 *   eventMs        —— MotionEvent 的 uptimeMillis（事件在系统输入子系统里"发生"的时刻）
 *   recvUptimeMs   —— 应用真正处理到它的时刻，同一 uptime 基准 → 两者之差就是**输入延迟**
 *   recvElapsedNs  —— 同一时刻的 elapsedRealtimeNanos，**与 SensorEvent.timestamp 同基准**
 *
 * 只有 recvElapsedNs 能和 IMU 的 CSV 直接对齐；输入延迟则是这次要顺带量出来的量。
 */
data class TouchMark(
    val name: String,
    val eventMs: Long,
    val recvUptimeMs: Long,
    val recvElapsedNs: Long,
)

/**
 * 落盘。格式沿用 docs/implementation-guide.md 阶段一：
 *
 * ```
 * acc.csv:  timestamp_ns,ax,ay,az
 * gyr.csv:  timestamp_ns,gx,gy,gz
 * meta.txt: 设备/传感器档案 + 体检报告 + 触摸事件
 * ```
 *
 * 两条流各自独立成文件，因为它们的回调时刻本来就不重合（时间对齐留给算法层）。
 * 写在 app 专属外部目录（`/sdcard/Android/data/<pkg>/files/slideruler/...`），
 * 不需要任何存储权限；取回用 `adb pull`（见 android/README.md）。
 */
object CaptureStore {

    fun save(
        context: Context,
        acc: List<Sample>,
        gyr: List<Sample>,
        accLabel: String,
        gyrLabel: String,
        report: String,
        marks: List<TouchMark> = emptyList(),
    ): File {
        val base = context.getExternalFilesDir(null) ?: context.filesDir
        val dir = File(File(base, "slideruler"), stamp()).apply { mkdirs() }
        writeCsv(File(dir, "acc.csv"), "timestamp_ns,ax,ay,az", acc)
        writeCsv(File(dir, "gyr.csv"), "timestamp_ns,gx,gy,gz", gyr)
        File(dir, "meta.txt").writeText(meta(acc, gyr, accLabel, gyrLabel, report, marks), Charsets.UTF_8)
        return dir
    }

    private fun writeCsv(file: File, header: String, samples: List<Sample>) {
        BufferedWriter(file.writer(Charsets.UTF_8), 1 shl 16).use { w ->
            w.append(header)
            w.newLine()
            val sb = StringBuilder(48)
            for (s in samples) {
                sb.setLength(0)
                sb.append(s.tNs).append(',')
                sb.append(fmt(s.x)).append(',')
                sb.append(fmt(s.y)).append(',')
                sb.append(fmt(s.z))
                w.append(sb)
                w.newLine()
            }
        }
    }

    private fun meta(
        acc: List<Sample>,
        gyr: List<Sample>,
        accLabel: String,
        gyrLabel: String,
        report: String,
        marks: List<TouchMark>,
    ): String = buildString {
        appendLine("# slideruler 采集档案")
        appendLine("时间      ${SimpleDateFormat("yyyy-MM-dd HH:mm:ss", Locale.US).format(Date())}")
        appendLine("设备      ${Build.MANUFACTURER} ${Build.MODEL}")
        appendLine("系统      Android ${Build.VERSION.RELEASE} (API ${Build.VERSION.SDK_INT})")
        appendLine("ABI       ${Build.SUPPORTED_ABIS.joinToString(", ")}")
        val requestedHz = 1_000_000.0 / ImuRecorder.DEFAULT_PERIOD_US
        appendLine("请求周期  ${ImuRecorder.DEFAULT_PERIOD_US} µs (${String.format(Locale.US, "%.0f", requestedHz)} Hz)")
        appendLine("加速度计  ${accLabel.replace('\n', ' ')}")
        appendLine("陀螺仪    ${gyrLabel.replace('\n', ' ')}")
        appendLine("样本数    acc ${acc.size} / gyr ${gyr.size}")
        appendLine("丢样阈值  dt > ${ImuRecorder.DROP_MS.toInt()} ms")
        if (marks.isNotEmpty()) {
            val t0 = acc.firstOrNull()?.tNs ?: 0L
            appendLine()
            appendLine("触摸事件（按住式边界）")
            appendLine("  recv_elapsed_ns 与 acc.csv 的时间戳同基准，可直接对齐；")
            appendLine("  输入延迟 = recv_uptime_ms - event_ms（同一 uptime 基准，仅 ms 精度）")
            for (m in marks) {
                val rel = (m.recvElapsedNs - t0) / 1e6
                appendLine(
                    String.format(
                        Locale.US, "  %-4s event=%d ms  recv_uptime=%d ms  recv_elapsed=%d ns  " +
                            "相对首样本=%.3f ms  输入延迟=%d ms",
                        m.name, m.eventMs, m.recvUptimeMs, m.recvElapsedNs, rel,
                        m.recvUptimeMs - m.eventMs,
                    )
                )
            }
        }
        appendLine()
        appendLine(report)
    }

    private fun stamp(): String =
        SimpleDateFormat("yyyyMMdd_HHmmss", Locale.US).format(Date())

    /** 不用 "%.6f".format：默认 Locale 可能把小数点写成逗号，破坏 CSV。 */
    private fun fmt(v: Float): String = String.format(Locale.US, "%.6f", v)
}
