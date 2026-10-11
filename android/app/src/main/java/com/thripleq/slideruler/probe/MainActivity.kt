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
 * 松手后再采多久才停。用来把"终点静止段"纳入缓冲区 ——
 * 按钮方案的终点边界不需要精确（落在静止段内即安全），但静止段本身要采到。
 */
private const val SLIDE_TAIL_MS = 700L

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
    val slideMarks: Int = 0,
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
            recording = true, slideMode = true, pressed = false, slideMarks = 0,
            mode = "滑动测量（按住式）", totalS = 0, remainS = 0,
            report = null, savedDir = null, error = null,
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
        ui = ui.copy(pressed = false, slideMarks = marks.size)
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
                recording = false, slideMode = false, pressed = false,
                error = "样本太少（acc ${acc.size} / gyr ${gyr.size}），本次不保存",
            )
            return
        }

        val report = buildSlideReport(acc, gyr.size)
        val dir = CaptureStore.save(
            this, acc, gyr, recorder.accLabel, recorder.gyrLabel, report, marks.toList(),
        )
        ui = ui.copy(
            recording = false, slideMode = false, pressed = false, error = null,
            report = report, savedDir = dir.absolutePath,
            accCount = acc.size, gyrCount = gyr.size,
        )
    }

    /** 用户取消本次滑动测量（未松手就退出 / 按了取消）。 */
    private fun exitSlide() {
        slidePending = false
        main.removeCallbacksAndMessages(null)
        if (recorder.isRecording) recorder.stop()
        window.clearFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON)
        ui = ui.copy(
            recording = false, slideMode = false, pressed = false,
            totalS = 0, remainS = 0, error = null,
        )
    }

    private fun buildSlideReport(acc: List<Sample>, gyrCount: Int): String {
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

            if (ui.slideMode) {
                Text(
                    "① 先静置 2 秒（学本机噪声底）→ ② 按住 → ③ 稍停 → ④ 推 → " +
                        "⑤ 停稳，等指示灯变绿 → ⑥ 立刻松手，别再碰手机。\n" +
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
                    if (ui.pressed) {
                        if (ui.still) "● 测量中 · 已静止，可以松手" else "● 测量中 · 还在动，别松手"
                    } else {
                        "按住开始"
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
