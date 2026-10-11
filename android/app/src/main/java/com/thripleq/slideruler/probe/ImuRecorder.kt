package com.thripleq.slideruler.probe

import android.content.Context
import android.hardware.Sensor
import android.hardware.SensorEvent
import android.hardware.SensorEventListener
import android.hardware.SensorManager
import kotlin.math.sqrt

/** 一条原始采样：纳秒时间戳 + 三轴原始值（加速度计 m/s²，陀螺 rad/s）。 */
data class Sample(val tNs: Long, val x: Float, val y: Float, val z: Float)

/**
 * 采集层。只做一件事：把两个**原始硬件传感器**的回调按各自的时间戳落进缓冲区。
 *
 * 刻意不用 TYPE_LINEAR_ACCELERATION / TYPE_GRAVITY —— 那是厂商黑盒派生量，质量不可控，
 * 官方文档明说可用性因设备而异。去重力我们自己做（见 docs/implementation-guide.md §1.2 第 6 条）。
 *
 * 两条流各自独立回调、采样时刻不重合（加速度计 12:00:00.010，陀螺 12:00:00.013）。
 * 这不是 bug 是现实，时间对齐留给算法层。
 */
class ImuRecorder(context: Context) {

    private val sm = context.getSystemService(Context.SENSOR_SERVICE) as SensorManager
    private val accSensor: Sensor? = sm.getDefaultSensor(Sensor.TYPE_ACCELEROMETER)
    private val gyrSensor: Sensor? = sm.getDefaultSensor(Sensor.TYPE_GYROSCOPE)

    private val accBuf = ArrayList<Sample>(65_536)
    private val gyrBuf = ArrayList<Sample>(65_536)

    private var lastAcc: Sample? = null
    private var lastGyr: Sample? = null
    private var listening = false

    val available: Boolean get() = accSensor != null && gyrSensor != null

    val accLabel: String get() = label(accSensor, "加速度计")
    val gyrLabel: String get() = label(gyrSensor, "陀螺仪")

    private fun label(s: Sensor?, fallback: String): String = when {
        s == null -> "无$fallback"
        else -> "${s.name}\n${s.vendor} · 最快 ${fastestHz(s.minDelay)} Hz · 量程 ${s.maximumRange}"
    }

    private fun fastestHz(minDelayUs: Int): String =
        if (minDelayUs > 0) "%.0f".format(1_000_000.0 / minDelayUs) else "?"

    private val listener = object : SensorEventListener {
        override fun onSensorChanged(event: SensorEvent) {
            val s = Sample(event.timestamp, event.values[0], event.values[1], event.values[2])
            when (event.sensor.type) {
                Sensor.TYPE_ACCELEROMETER -> {
                    if (accBuf.size < MAX_SAMPLES) accBuf.add(s)
                    lastAcc = s
                }
                Sensor.TYPE_GYROSCOPE -> {
                    if (gyrBuf.size < MAX_SAMPLES) gyrBuf.add(s)
                    lastGyr = s
                }
            }
        }

        override fun onAccuracyChanged(sensor: Sensor?, accuracy: Int) = Unit
    }

    /** 开始采集：清空缓冲，按 [periodUs] 请求两个传感器。返回是否成功。 */
    fun start(periodUs: Int = DEFAULT_PERIOD_US): Boolean {
        if (!available) return false
        accBuf.clear()
        gyrBuf.clear()
        lastAcc = null
        lastGyr = null
        sm.registerListener(listener, accSensor, periodUs)
        sm.registerListener(listener, gyrSensor, periodUs)
        listening = true
        return true
    }

    fun stop() {
        if (listening) sm.unregisterListener(listener)
        listening = false
    }

    val isRecording: Boolean get() = listening

    fun accSamples(): List<Sample> = accBuf
    fun gyrSamples(): List<Sample> = gyrBuf

    val accLast: Sample? get() = lastAcc
    val gyrLast: Sample? get() = lastGyr

    /** 最近一次相邻采样的间隔（ms）——实时显示 jitter 用。 */
    fun lastAccDtMs(): Double = tailDtMs(accBuf)
    fun lastGyrDtMs(): Double = tailDtMs(gyrBuf)

    /**
     * 最近 [win] 个样本的「运动强度」，给界面做一个实时稳定度指示：
     *
     *   accStd = 三轴标准差的较大者 —— 不需要姿态解算，最鲁棒（指南 §7 判据三）
     *   gyrMax = 陀螺范数峰值（指南 §7 判据二）
     *
     * 静止时两者都只反映本机噪声地板（K40 实测 σ_a≈0.019、σ_g≈0.006）；
     * 一旦被手碰到或在移动，两者会同时抬头。**不用「‖f‖ 偏离 9.81」那条判据**——
     * 仿真已证明它横向几乎瞎（见 sim/README.md 勘误 ①）。
     */
    fun motionStats(win: Int = 40): Pair<Double, Double> {
        val a = accBuf
        val g = gyrBuf
        if (a.size < 5 || g.size < 5) return 0.0 to 0.0
        val na = if (win < a.size) win else a.size
        var sx = 0.0
        var sy = 0.0
        var sz = 0.0
        for (k in a.size - na until a.size) {
            sx += a[k].x.toDouble()
            sy += a[k].y.toDouble()
            sz += a[k].z.toDouble()
        }
        sx /= na
        sy /= na
        sz /= na
        var vx = 0.0
        var vy = 0.0
        var vz = 0.0
        for (k in a.size - na until a.size) {
            val dx = a[k].x - sx
            val dy = a[k].y - sy
            val dz = a[k].z - sz
            vx += dx * dx
            vy += dy * dy
            vz += dz * dz
        }
        val denom = (na - 1).coerceAtLeast(1)
        var best = vx
        if (vy > best) best = vy
        if (vz > best) best = vz
        val accStd = sqrt(best / denom)

        val ng = if (win < g.size) win else g.size
        var gmax = 0.0
        for (k in g.size - ng until g.size) {
            val x = g[k].x.toDouble()
            val y = g[k].y.toDouble()
            val z = g[k].z.toDouble()
            val m = sqrt(x * x + y * y + z * z)
            if (m > gmax) gmax = m
        }
        return accStd to gmax
    }

    private fun tailDtMs(b: List<Sample>): Double =
        if (b.size >= 2) (b[b.size - 1].tNs - b[b.size - 2].tNs) / 1e6 else 0.0

    companion object {
        /** 10 000 µs = 100 Hz，与 docs/implementation-guide.md 阶段一一致。 */
        const val DEFAULT_PERIOD_US = 10_000

        /** 采样间隔超过该值视为一次丢样/调度断续。 */
        const val DROP_MS = 50.0

        /** 100 Hz 下约 50 分钟，纯粹防 OOM。 */
        const val MAX_SAMPLES = 300_000
    }
}
