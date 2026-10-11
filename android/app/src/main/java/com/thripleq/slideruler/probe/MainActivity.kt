package com.thripleq.slideruler.probe

import android.os.Build
import android.os.Bundle
import android.os.Handler
import android.os.Looper
import android.os.SystemClock
import android.view.WindowManager
import androidx.activity.ComponentActivity
import androidx.activity.compose.setContent
import androidx.activity.enableEdgeToEdge
import androidx.compose.foundation.background
import androidx.compose.foundation.gestures.awaitEachGesture
import androidx.compose.foundation.gestures.awaitFirstDown
import androidx.compose.foundation.layout.Arrangement
import androidx.compose.foundation.layout.Box
import androidx.compose.foundation.layout.Column
import androidx.compose.foundation.layout.Row
import androidx.compose.foundation.layout.fillMaxSize
import androidx.compose.foundation.layout.fillMaxWidth
import androidx.compose.foundation.layout.height
import androidx.compose.foundation.layout.padding
import androidx.compose.foundation.layout.size
import androidx.compose.foundation.rememberScrollState
import androidx.compose.foundation.shape.CircleShape
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.foundation.verticalScroll
import androidx.compose.material3.Button
import androidx.compose.material3.Card
import androidx.compose.material3.LinearProgressIndicator
import androidx.compose.material3.MaterialTheme
import androidx.compose.material3.OutlinedButton
import androidx.compose.material3.Scaffold
import androidx.compose.material3.Surface
import androidx.compose.material3.Text
import androidx.compose.runtime.Composable
import androidx.compose.runtime.getValue
import androidx.compose.runtime.mutableStateOf
import androidx.compose.runtime.setValue
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.input.pointer.pointerInput
import androidx.compose.ui.text.font.FontFamily
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import java.util.Locale
import kotlin.math.sqrt

/** 体检 A 的静止时长，对应 docs/implementation-guide.md §2。 */
private const val HEALTH_SECONDS = 30
private const val TICK_MS = 100L

/**
 * 松手后再采多久才停。
 *
 * **必须 ≥ 引擎的停稳搜索上限（SlideMeasure.MAX_EXT_S = 1.2 s）**：松手时手机还在动的话，
 * 终点是算法在松手之后才找到的，缓冲区里得先有那段数据。留 0.3 s 余量 ——
 * 换算成用户看到的规则就是：松手后保持手机不动 1.5 s 才出数。
 */
private const val SLIDE_TAIL_MS = 1500L

/**
 * 稳定度指示的基线学习时长。取 2.5 s，和引导文案「先静置 2 秒」对齐：
 * 前 2.5 s 取到的最小运动量当作这台机当次的噪声地板，之后冻结。
 */
private const val STAB_CAL_MS = 2500L

/**
 * 静止判据的绝对下限（跨设备兜底，防止"基线学到运动量"导致永远显示静止）。
 * acc 标准差下限取 K40 实测噪声 σ_a≈0.019 的 ~5 倍；陀螺范数下限同理。
 */
private const val STAB_ACC_FLOOR = 0.10
private const val STAB_GYR_FLOOR = 0.03

data class UiState(
    val device: String = "",
    val accLabel: String = "",
    val gyrLabel: String = "",
    val available: Boolean = true,
    val recording: Boolean = false,
    val mode: String = "",
    val totalS: Int = 0,
    val remainS: Int = 0,
    val acc: Sample? = null,
    val gyr: Sample? = null,
    val accCount: Int = 0,
    val gyrCount: Int = 0,
    val accRateHz: Double = 0.0,
    val gyrRateHz: Double = 0.0,
    val accDtMs: Double = 0.0,
    val report: String? = null,
    val savedDir: String? = null,
    val error: String? = null,
    // ---- 滑动测量（按住式边界）----
    val slideMode: Boolean = false,
    val pressed: Boolean = false,
    /** 已松手、正在补采尾部静止段。这 1.5 s 里别动手机，算法要靠它定终点。 */
    val slideTail: Boolean = false,
    val slideMarks: Int = 0,
    /** 本次测出的距离（cm）。null = 还没测或没测出来。 */
    val slideDistanceCm: Double? = null,
    /** 一句话结论：这一次可不可信。 */
    val slideVerdict: String? = null,
    /** 距离下面的诊断细节（窗口 / 判据 / 本机参数）。 */
    val slideDetail: String? = null,
    // ---- 实时稳定度指示 ----
    val still: Boolean = false,
    val stabScore: Double = 0.0,
    val stabAcc: Double = 0.0,
    val stabGyr: Double = 0.0,
)

class MainActivity : ComponentActivity() {

    private lateinit var recorder: ImuRecorder
    private var ui by mutableStateOf(UiState())

    private val main = Handler(Looper.getMainLooper())
    private var startMs = 0L
    private var durationS = 0

    /** 滑动测量的触摸标记（DOWN/UP），落盘时写进 meta.txt。 */
    private val marks = ArrayList<TouchMark>()
    /** 已收到 UP、正在等尾部静止段采完，期间忽略新的触摸。 */
    private var slidePending = false
    /** 上一次测距引擎的耗时（ms）——纯诊断，放进结果里给用户看。 */
    private var lastMeasureMs = 0L

    /**
     * 稳定度指示用的「本机噪声地板」——录制过程中持续取到的最小值。
     * 手机真的静止时，[ImuRecorder.motionStats] 会掉到地板；被手碰到就抬头。
     */
    private var stabMinAcc = Double.MAX_VALUE
    private var stabMinGyr = Double.MAX_VALUE
    /** 基线只在前 [STAB_CAL_MS] 学，之后冻结 —— 否则运行最小值会一直下探、阈值越来越紧。 */
    private var stabFrozen = false

    private val ticker = object : Runnable {
        override fun run() {
            val elapsedMs = SystemClock.elapsedRealtime() - startMs
            val elapsedS = elapsedMs / 1000.0
            val remain = if (durationS > 0) (durationS - elapsedS).toInt().coerceAtLeast(0) else 0
            val ac = recorder.accSamples().size
            val gc = recorder.gyrSamples().size

            // 稳定度：阈值 = max(4 × 本机地板, 绝对下限)，取两者更差的那个比值
            val (accStd, gyrMax) = recorder.motionStats()
            if (accStd > 0.0 && !stabFrozen) {
                if (accStd < stabMinAcc) stabMinAcc = accStd
                if (gyrMax < stabMinGyr) stabMinGyr = gyrMax
                if (elapsedMs >= STAB_CAL_MS) stabFrozen = true
            }
            val thrA = maxOf(stabMinAcc * 4.0, STAB_ACC_FLOOR)
            val thrG = maxOf(stabMinGyr * 4.0, STAB_GYR_FLOOR)
            val score = maxOf(accStd / thrA, gyrMax / thrG)

            ui = ui.copy(
                acc = recorder.accLast,
                gyr = recorder.gyrLast,
                accCount = ac,
                gyrCount = gc,
                accDtMs = recorder.lastAccDtMs(),
                accRateHz = if (elapsedS > 1.0) ac / elapsedS else 0.0,
                gyrRateHz = if (elapsedS > 1.0) gc / elapsedS else 0.0,
                remainS = remain,
                still = score < 1.0,
                stabScore = score,
                stabAcc = accStd,
                stabGyr = gyrMax,
            )
            if (durationS > 0 && remain <= 0) {
                finishCapture()
                return
            }
            if (ui.recording) main.postDelayed(this, TICK_MS)
        }
    }

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        enableEdgeToEdge()
        recorder = ImuRecorder(this)
        ui = ui.copy(
            device = "${Build.MANUFACTURER} ${Build.MODEL} · Android ${Build.VERSION.RELEASE} (API ${Build.VERSION.SDK_INT})",
            accLabel = recorder.accLabel,
            gyrLabel = recorder.gyrLabel,
            available = recorder.available,
        )
        setContent {
            MaterialTheme {
                Surface(Modifier.fillMaxSize(), color = MaterialTheme.colorScheme.background) {
                    ProbeScreen(
                        ui = ui,
                        onHealth = { startCapture(HEALTH_SECONDS, "体检 A（静止 ${HEALTH_SECONDS} s）") },
                        onFree = { startCapture(0, "自由采集") },
                        onStop = { finishCapture() },
                        onSlideEnter = { enterSlide() },
                        onSlideDown = { e, u, n -> slideDown(e, u, n) },
                        onSlideUp = { e, u, n -> slideUp(e, u, n) },
                        onSlideExit = { exitSlide() },
                    )
                }
            }
        }
    }

    private fun startCapture(seconds: Int, mode: String) {
        if (!recorder.available) {
            ui = ui.copy(error = "本机缺少加速度计或陀螺仪，无法采集")
            return
        }
        if (!recorder.start()) {
            ui = ui.copy(error = "传感器注册失败")
            return
        }
        stabMinAcc = Double.MAX_VALUE
        stabMinGyr = Double.MAX_VALUE
        stabFrozen = false
        durationS = seconds
        startMs = SystemClock.elapsedRealtime()
        window.addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON)
        ui = ui.copy(
            recording = true, mode = mode, totalS = seconds, remainS = seconds,
            report = null, savedDir = null, error = null,
        )
        main.removeCallbacks(ticker)
        main.postDelayed(ticker, TICK_MS)
    }

    private fun finishCapture() {
        main.removeCallbacks(ticker)
        val wasRecording = ui.recording || recorder.isRecording
        recorder.stop()
        window.clearFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON)
        if (!wasRecording) return

        val acc = recorder.accSamples().toList()
        val gyr = recorder.gyrSamples().toList()
        if (acc.size < 100 || gyr.size < 100) {
            ui = ui.copy(
                recording = false, totalS = 0, remainS = 0,
                error = "样本太少（acc ${acc.size} / gyr ${gyr.size}），本次不保存",
            )
            return
        }

        val report = HealthCheck.buildReport(
            HealthCheck.analyse(acc),
            HealthCheck.analyse(gyr),
            ui.mode,
        )
        val dir = CaptureStore.save(this, acc, gyr, recorder.accLabel, recorder.gyrLabel, report)
        ui = ui.copy(
            recording = false, totalS = 0, remainS = 0, error = null,
            report = report, savedDir = dir.absolutePath,
            accCount = acc.size, gyrCount = gyr.size,
        )
    }

    // ---------------------------------------------------------------- 滑动测量

    /** 进入滑动测量：立刻开始连续采集，等用户按住按钮来划定边界。 */
    private fun enterSlide() {
        if (!recorder.available) {
            ui = ui.copy(error = "本机缺少加速度计或陀螺仪，无法采集")
            return
        }
        if (!recorder.start()) {
            ui = ui.copy(error = "传感器注册失败")
            return
        }
        marks.clear()
        slidePending = false
        stabMinAcc = Double.MAX_VALUE
        stabMinGyr = Double.MAX_VALUE
        stabFrozen = false
        durationS = 0
        startMs = SystemClock.elapsedRealtime()
        window.addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON)
        ui = ui.copy(
            recording = true, slideMode = true, pressed = false, slideTail = false, slideMarks = 0,
            mode = "滑动测量（按住式）", totalS = 0, remainS = 0,
            report = null, savedDir = null, error = null,
            slideDistanceCm = null, slideVerdict = null, slideDetail = null,
        )
        main.removeCallbacks(ticker)
        main.postDelayed(ticker, TICK_MS)
    }

    private fun slideDown(eventMs: Long, recvUptimeMs: Long, recvElapsedNs: Long) {
        if (!ui.slideMode || slidePending) return
        if (marks.isNotEmpty()) return // 本次已有标记（等保存或已在测量中）
        marks.add(TouchMark("DOWN", eventMs, recvUptimeMs, recvElapsedNs))
        ui = ui.copy(pressed = true, slideMarks = marks.size)
    }

    private fun slideUp(eventMs: Long, recvUptimeMs: Long, recvElapsedNs: Long) {
        if (!ui.slideMode || slidePending) return
        if (marks.none { it.name == "DOWN" }) return
        if (marks.any { it.name == "UP" }) return
        marks.add(TouchMark("UP", eventMs, recvUptimeMs, recvElapsedNs))
        slidePending = true
        ui = ui.copy(pressed = false, slideTail = true, slideMarks = marks.size)
        // 再采一小段静止，让终点边界之外也有静止段
        main.postDelayed({ finishSlide() }, SLIDE_TAIL_MS)
    }

    private fun finishSlide() {
        slidePending = false
        main.removeCallbacks(ticker)
        val wasRecording = recorder.isRecording
        recorder.stop()
        window.clearFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON)
        if (!wasRecording) return

        val acc = recorder.accSamples().toList()
        val gyr = recorder.gyrSamples().toList()
        if (acc.size < 100 || gyr.size < 100) {
            ui = ui.copy(
                recording = false, slideMode = false, pressed = false, slideTail = false,
                error = "样本太少（acc ${acc.size} / gyr ${gyr.size}），本次不保存",
            )
            return
        }

        // 设备端出数。以前距离只能 adb pull 回 PC 再算 —— 那就不是一把尺子。
        val result = measureSlide(acc, gyr)
        val report = buildSlideReport(acc, gyr.size, result)
        val dir = CaptureStore.save(
            this, acc, gyr, recorder.accLabel, recorder.gyrLabel, report, marks.toList(),
        )
        ui = ui.copy(
            recording = false, slideMode = false, pressed = false, slideTail = false,
            error = result?.takeIf { !it.ok }?.let { "测距失败：${it.failReason}" },
            report = report, savedDir = dir.absolutePath,
            accCount = acc.size, gyrCount = gyr.size,
            slideDistanceCm = result?.takeIf { it.ok }?.distanceCm,
            slideVerdict = result?.let { slideVerdict(it) },
            slideDetail = result?.let { slideDetail(it) },
        )
    }

    /** 跑一次测距引擎。没拿到完整的 DOWN/UP 就不测（返回 null）。 */
    private fun measureSlide(acc: List<Sample>, gyr: List<Sample>): MeasureResult? {
        val down = marks.firstOrNull { it.name == "DOWN" } ?: return null
        val up = marks.firstOrNull { it.name == "UP" } ?: return null
        val t0 = SystemClock.elapsedRealtime()
        val r = SlideMeasure.measure(
            acc.map { it.toPt3() },
            gyr.map { it.toPt3() },
            down.recvElapsedNs,
            up.recvElapsedNs,
        )
        lastMeasureMs = SystemClock.elapsedRealtime() - t0
        return r
    }

    /** 一句话结论 —— 用户需要知道这一次是怎么来的。 */
    private fun slideVerdict(r: MeasureResult): String {
        if (!r.ok) return r.failReason ?: "测距失败"
        return if (r.extended) {
            "松手时手机还在动 → 算法往前找到真正停稳的位置，补回 " +
                String.format(Locale.US, "%.0f", r.extendedMs) + " ms"
        } else {
            "松手时手机已经停稳（终点速度在阈值内，无需顺延）"
        }
    }

    /** 距离下面的诊断：窗口、判据、本机参数。用户看着它判断"这次为什么是这个数"。 */
    private fun slideDetail(r: MeasureResult): String {
        if (!r.ok) return "原因：${r.failReason}"
        val rel = if (kotlin.math.abs(r.vEndFwd) <= r.threshold) "≤" else ">"
        return buildString {
            appendLine(
                String.format(
                    Locale.US, "窗口 #%d→#%d（按住 %.2f s / 积分 %.2f s / 顺延 %.0f ms）",
                    r.downIndex, r.endIndex, r.holdS, r.durationS, r.extendedMs,
                )
            )
            appendLine(
                String.format(
                    Locale.US, "终点速度 %+.4f m/s %s 阈值 %.4f → %s",
                    r.vEndFwd, rel, r.threshold, if (r.detrended) "矢量去趋势" else "裸积分",
                )
            )
            appendLine(
                String.format(
                    Locale.US, "本机 |g| %.4f   σ_a %.4f   |b| %.4f   dt %.3f ms",
                    r.gravityMag, r.accSigma, r.biasBound, r.dt * 1000.0,
                )
            )
            append(
                String.format(
                    Locale.US, "裸积分对照 %.2f cm   引擎耗时 %d ms",
                    r.rawDistanceCm, lastMeasureMs,
                )
            )
        }
    }

    /** 用户取消本次滑动测量（未松手就退出 / 按了取消）。 */
    private fun exitSlide() {
        slidePending = false
        main.removeCallbacksAndMessages(null)
        if (recorder.isRecording) recorder.stop()
        window.clearFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON)
        ui = ui.copy(
            recording = false, slideMode = false, pressed = false, slideTail = false,
            totalS = 0, remainS = 0, error = null,
        )
    }

    private fun buildSlideReport(acc: List<Sample>, gyrCount: Int, result: MeasureResult?): String {
        if (acc.isEmpty()) return "（无样本）"
        val t0 = acc.first().tNs
        val sb = StringBuilder()
        sb.appendLine("滑动测量（按住式边界）")
        sb.appendLine("样本  acc ${acc.size} / gyr $gyrCount")
        if (marks.isEmpty()) {
            sb.appendLine("未捕获到触摸事件")
            return sb.toString()
        }
        val down = marks.firstOrNull { it.name == "DOWN" }
        val up = marks.firstOrNull { it.name == "UP" }
        for (m in marks) {
            val rel = (m.recvElapsedNs - t0) / 1e6
            sb.appendLine(
                String.format(
                    Locale.US, "%s  相对首样本 %.1f ms  输入延迟 %d ms",
                    m.name, rel, m.recvUptimeMs - m.eventMs,
                )
            )
        }
        if (down != null && up != null) {
            val spanMs = (up.eventMs - down.eventMs)
            sb.appendLine(String.format(Locale.US, "按住时长（事件时刻之差） %d ms", spanMs))
        }
        when {
            result == null -> sb.appendLine("未捕获到完整的 DOWN/UP —— 本次无法测距")
            !result.ok -> sb.appendLine("测距失败：${result.failReason}")
            else -> {
                sb.appendLine(
                    String.format(
                        Locale.US, "★ 距离 %.2f cm（裸积分 %.2f cm）",
                        result.distanceCm, result.rawDistanceCm,
                    )
                )
                sb.appendLine("  " + slideVerdict(result))
                sb.appendLine("  " + slideDetail(result).replace('\n', ' ').trim())
            }
        }
        return sb.toString()
    }

    override fun onDestroy() {
        main.removeCallbacksAndMessages(null)
        if (::recorder.isInitialized) recorder.stop()
        super.onDestroy()
    }
}

// ------------------------------------------------------------------ 界面

@Composable
private fun ProbeScreen(
    ui: UiState,
    onHealth: () -> Unit,
    onFree: () -> Unit,
    onStop: () -> Unit,
    onSlideEnter: () -> Unit,
    onSlideDown: (Long, Long, Long) -> Unit,
    onSlideUp: (Long, Long, Long) -> Unit,
    onSlideExit: () -> Unit,
) {
    Scaffold { pad ->
        Column(
            modifier = Modifier
                .fillMaxSize()
                .padding(pad)
                .verticalScroll(rememberScrollState())
                .padding(16.dp),
            verticalArrangement = Arrangement.spacedBy(12.dp),
        ) {
            Text("slideruler · 采集 / 体检 A", style = MaterialTheme.typography.titleMedium)

            Card(Modifier.fillMaxWidth()) {
                Column(Modifier.padding(12.dp), verticalArrangement = Arrangement.spacedBy(4.dp)) {
                    Text(ui.device, fontSize = 13.sp)
                    if (!ui.available) {
                        Text(
                            "✗ 本机缺少加速度计或陀螺仪",
                            color = MaterialTheme.colorScheme.error,
                            fontSize = 13.sp,
                        )
                    }
                    Mono(ui.accLabel, 11)
                    Mono(ui.gyrLabel, 11)
                }
            }

            Card(Modifier.fillMaxWidth()) {
                Column(Modifier.padding(12.dp), verticalArrangement = Arrangement.spacedBy(4.dp)) {
                    Text(
                        if (ui.recording) "采集中 · ${ui.mode}" else "空闲",
                        fontSize = 13.sp,
                    )
                    Mono("acc  ${vec(ui.acc)}   |a| ${mag(ui.acc)}", 12)
                    Mono("gyr  ${vec(ui.gyr)}   |ω| ${mag(ui.gyr)}", 12)
                    Mono(
                        "实测 acc ${hz(ui.accRateHz)} · gyr ${hz(ui.gyrRateHz)}" +
                            "   dt ${f(ui.accDtMs, 2)} ms",
                        12,
                    )
                    Mono("样本 acc ${ui.accCount} · gyr ${ui.gyrCount}", 12)
                }
            }

            val distanceCm = ui.slideDistanceCm
            if (distanceCm != null) {
                Card(Modifier.fillMaxWidth()) {
                    Column(
                        Modifier.padding(12.dp),
                        verticalArrangement = Arrangement.spacedBy(6.dp),
                    ) {
                        Text("距离", style = MaterialTheme.typography.titleSmall)
                        Text(
                            "${f(distanceCm, 2)} cm",
                            fontSize = 40.sp,
                            fontFamily = FontFamily.Monospace,
                            textAlign = TextAlign.Center,
                            modifier = Modifier.fillMaxWidth(),
                        )
                        ui.slideVerdict?.let { Text(it, fontSize = 13.sp) }
                        ui.slideDetail?.let {
                            Text(it, fontSize = 11.sp, fontFamily = FontFamily.Monospace)
                        }
                    }
                }
            }

            if (ui.slideMode) {
                Text(
                    "① 先静置 2 秒（学本机噪声底）→ ② 按住 → ③ 推 → ④ 松手 → " +
                        "⑤ 保持手机不动 1.5 秒 → 出结果。\n" +
                        "松手时还在动也不要紧：算法会往前找到手机真正停住的位置补回来。" +
                        "但松手后手机必须在 1.5 秒内停下 —— 提前把它拿起来就会算不回来。\n" +
                        "提示：按住到松手尽量控制在 1~2 秒内 —— 拖得越久误差越大。",
                    fontSize = 13.sp,
                )
                StabilityBar(ui)
                SlidePad(
                    pressed = ui.pressed,
                    onDown = onSlideDown,
                    onUp = onSlideUp,
                )
                Text(
                    when {
                        ui.slideTail -> "已松手 · 保持手机不动，正在算…"
                        ui.pressed && ui.still -> "● 测量中 · 已静止，可以松手"
                        ui.pressed -> "● 测量中 · 还在动，别松手"
                        else -> "按住开始"
                    },
                    fontSize = 13.sp,
                    textAlign = TextAlign.Center,
                    modifier = Modifier.fillMaxWidth(),
                )
                OutlinedButton(onClick = onSlideExit) { Text("取消本次") }
            } else {
                Row(horizontalArrangement = Arrangement.spacedBy(8.dp)) {
                    Button(onClick = onHealth, enabled = !ui.recording && ui.available) {
                        Text("体检 A · 静止 ${HEALTH_SECONDS} s")
                    }
                    if (ui.recording) {
                        OutlinedButton(onClick = onStop) { Text("停止并保存") }
                    } else {
                        Button(onClick = onFree, enabled = ui.available) { Text("自由采集") }
                    }
                }
                Button(
                    onClick = onSlideEnter,
                    enabled = !ui.recording && ui.available,
                    modifier = Modifier.fillMaxWidth(),
                ) {
                    Text("滑动测量（按住式边界）")
                }
            }

            if (ui.recording && !ui.slideMode && ui.totalS > 0) {
                val frac = ((ui.totalS - ui.remainS).toFloat() / ui.totalS).coerceIn(0f, 1f)
                LinearProgressIndicator(progress = { frac }, modifier = Modifier.fillMaxWidth())
                Text("剩余 ${ui.remainS} s", fontSize = 13.sp)
            } else if (ui.recording && !ui.slideMode) {
                Text("自由采集中 … 点「停止并保存」结束", fontSize = 13.sp)
            }

            ui.error?.let {
                Text(it, color = MaterialTheme.colorScheme.error, fontSize = 13.sp)
            }

            val savedDir = ui.savedDir
            if (savedDir != null) {
                Text("已保存到：", fontSize = 12.sp)
                Mono(savedDir, 11)
            }

            ui.report?.let { report ->
                Card(Modifier.fillMaxWidth()) {
                    Column(Modifier.padding(12.dp)) {
                        Text("结果", style = MaterialTheme.typography.titleSmall)
                        Text(
                            report,
                            fontSize = 11.sp,
                            fontFamily = FontFamily.Monospace,
                            modifier = Modifier.padding(top = 6.dp),
                        )
                    }
                }
            }
        }
    }
}

/**
 * 实时稳定度指示。
 *
 * 判据：最近 0.4 s 的三轴加速度标准差 / 陀螺范数峰值，与「本机噪声地板 × 3」比较。
 * 立这条的理由是硬的 —— 真机第一次滑动（20261011_084700）失败的原因就是
 * 「松手时手机还在动」：算出来的距离随积分端点从 0 滑到 45 cm，没有平台。
 * 那个失效模式用户当场就能看见，所以必须把它显示出来。
 */
@Composable
private fun StabilityBar(ui: UiState) {
    val green = Color(0xFF2E7D32)
    val red = Color(0xFFC62828)
    Row(
        modifier = Modifier.fillMaxWidth(),
        verticalAlignment = Alignment.CenterVertically,
        horizontalArrangement = Arrangement.spacedBy(8.dp),
    ) {
        Box(
            Modifier
                .size(16.dp)
                .clip(CircleShape)
                .background(if (ui.still) green else red)
        )
        Text(
            if (ui.still) "静止 · 可以松手" else "在动 · 继续按住",
            fontSize = 13.sp,
            color = if (ui.still) green else red,
        )
        Mono(
            "aσ ${f(ui.stabAcc, 3)}  ω ${f(ui.stabGyr, 3)}  负荷 ${f(ui.stabScore, 2)}×",
            10,
        )
    }
}

/**
 * 按住式测量的大按钮。
 *
 * 直接用底层指针事件而不是 `clickable` / `detectTapGestures`，因为要拿到
 * `PointerInputChange.uptimeMillis`（事件**发生**的时刻），而不仅仅是"回调被调用的时刻"。
 * 两者之差就是输入延迟 —— 这次要顺带量出来的量。
 */
@Composable
private fun SlidePad(
    pressed: Boolean,
    onDown: (Long, Long, Long) -> Unit,
    onUp: (Long, Long, Long) -> Unit,
) {
    val idle = MaterialTheme.colorScheme.surfaceVariant
    val active = MaterialTheme.colorScheme.primary
    Box(
        modifier = Modifier
            .fillMaxWidth()
            .height(180.dp)
            .clip(RoundedCornerShape(18.dp))
            .background(if (pressed) active else idle)
            .pointerInput(Unit) {
                awaitEachGesture {
                    val down = awaitFirstDown(requireUnconsumed = false)
                    onDown(
                        down.uptimeMillis,
                        SystemClock.uptimeMillis(),
                        SystemClock.elapsedRealtimeNanos(),
                    )
                    while (true) {
                        val ev = awaitPointerEvent()
                        val ch = ev.changes.firstOrNull { it.id == down.id } ?: break
                        if (!ch.pressed) {
                            onUp(
                                ch.uptimeMillis,
                                SystemClock.uptimeMillis(),
                                SystemClock.elapsedRealtimeNanos(),
                            )
                            break
                        }
                    }
                }
            },
        contentAlignment = Alignment.Center,
    ) {
        Text(
            if (pressed) "● 测量中" else "按住开始",
            fontSize = 22.sp,
            color = if (pressed) MaterialTheme.colorScheme.onPrimary
            else MaterialTheme.colorScheme.onSurfaceVariant,
        )
    }
}

@Composable
private fun Mono(text: String, size: Int) {
    Text(text, fontSize = size.sp, fontFamily = FontFamily.Monospace)
}

private fun vec(s: Sample?): String =
    if (s == null) "--" else "${f(s.x.toDouble(), 3)}  ${f(s.y.toDouble(), 3)}  ${f(s.z.toDouble(), 3)}"

private fun mag(s: Sample?): String {
    if (s == null) return "--"
    val x = s.x.toDouble()
    val y = s.y.toDouble()
    val z = s.z.toDouble()
    return f(sqrt(x * x + y * y + z * z), 3)
}

private fun hz(v: Double): String = if (v <= 0.0) "--" else "${f(v, 1)} Hz"

/** 固定用 US Locale：默认 Locale 可能把小数点写成逗号。 */
private fun f(v: Double, digits: Int): String = String.format(Locale.US, "%.${digits}f", v)
