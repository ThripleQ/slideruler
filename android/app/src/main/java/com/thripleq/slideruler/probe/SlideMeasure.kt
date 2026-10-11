package com.thripleq.slideruler.probe

import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt

/**
 * 测距引擎：把 `sim/slide_pipeline.py` 的链路原样搬到设备端。
 *
 * 设备端本来只负责采集，距离要靠 `adb pull` 回 PC 再算 —— 那就不是一把尺子。
 * 这个文件补上最后一环：松手后直接出数。它是纯逻辑（不碰任何 Android API），
 * 所以能用 JVM 单测拿真机采集逐位对齐 Python 的结果（见 SlideMeasureTest）。
 *
 * 链路（与 Python 一一对应）：
 *
 *   1. 时间基准   dt 取「相邻 acc 时间戳间隔的中位数」—— 不是 10 ms。
 *                 位移 ∝ dt²，写死 10 ms 会在 K40 上凭空多出 11%（见 docs/calibration.md）。
 *   2. 边界       DOWN/UP 的 recv_elapsed_ns 与 Sample.tNs 同基准，二分落到样本下标。
 *   3. 自校准     只用 DOWN 之前 2 s 的静止段，估陀螺零偏 / 重力参考 / 噪声地板。
 *                 **不把任何常数写死**，每台机每次会话重算。
 *   4. 姿态       四元数互补滤波（α=0.05，振动超 0.3 rad 时只外推不修正）。
 *   5. 去重力     扣的是**本机实测的重力模长**，不是 9.81。
 *                 K40 实测 ‖g‖=9.9129，与 9.81 差 0.10 m/s² —— 等于噪声地板的 12 倍，
 *                 慢速滑动时会直接把主轴带偏到重力方向（详见 docs/boundary-design.md §3.8）。
 *   6. 积分       三维梯形二次积分，取**净位移矢量的模**。
 *                 不用 PCA 主轴 + 一维投影：投影方向取的是"加速度方差最大"，
 *                 而我们要的是"净位移方向"，真机上两者差 16°~24°，自带 cos 的损失。
 *   7. 条件 ZUPT  |v_end| 大到零偏解释不了时才去趋势；去趋势按三轴各自扣线性漂移。
 *   8. 终点顺延   松手时手机还在动 → 往后找它真正停稳的那一瞬，把终点挪过去。
 *                 （判据与第 7 步共用同一个阈值，两个决定永远自洽。）
 *
 * 刻意不做的两件事：
 *   · 不开「单姿态标度修正」k̂=‖g‖/9.81 —— 单姿态下 k̂ 混着零偏在重力方向的投影
 *     b·ĝ，分不开（b_z 只要 0.05 m/s² 就装成 +0.5% 的标度）。见 PipelineParams.use_scale_corr。
 *   · 不做静止检测器自动定位窗口 —— 按钮边界已经够准，检测器只当对照。
 */

// ------------------------------------------------------------------ 参数

private const val ALPHA = 0.05           // 互补滤波系数
private const val REASON_OK_RAD = 0.3    // 振动防护：误差角速度超此值时只外推不修正
private const val CAL_MAX_S = 2.0        // 校准窗上限（DOWN 之前的静止段）
private const val K_SIGMA = 5.0          // 条件 ZUPT 阈值系数，thr = |b|·T + K·σ·√(T·dt)
private const val MAX_EXT_S = 1.2        // 松手时还在动 → 最多往后找 1.2 s 的停稳点
private const val DWELL_S = 0.05         // 「停稳」需要持续多久才算数
private const val MIN_SAMPLES = 40

/** 一条三轴样本。用 Double 而不是 Float，是为了让单测能与 Python 逐位对齐。 */
class Pt3(val tNs: Long, val x: Double, val y: Double, val z: Double)

/** 采集层的 [Sample]（Float）→ 引擎的 [Pt3]（Double）。 */
fun Sample.toPt3(): Pt3 = Pt3(tNs, x.toDouble(), y.toDouble(), z.toDouble())

/**
 * 测量结果。[ok] 为 false 时 [failReason] 给出原因，其余字段无意义。
 *
 * 除了距离，刻意把「引擎看到的证据」一并带出来 —— 用户需要知道这一次可不可信：
 * [extended] 表示松手时手机还在动、终点是算法自己找的；[detrended] 表示做了去趋势。
 */
class MeasureResult(
    val ok: Boolean,
    val failReason: String? = null,
    /** 净位移（米）。 */
    val distanceM: Double = 0.0,
    /** 三维裸积分（未去趋势），做对照用。 */
    val rawDistanceM: Double = 0.0,
    /** 沿前进方向的终点速度（m/s）。接近 0 才是"真的停了"。 */
    val vEndFwd: Double = 0.0,
    /** 松手判据的阈值（m/s）。 */
    val threshold: Double = 0.0,
    val detrended: Boolean = false,
    /** 是否发生了终点顺延（松手时手机还在动）。 */
    val extended: Boolean = false,
    val extendedMs: Double = 0.0,
    val downIndex: Int = 0,
    val upIndex: Int = 0,
    val endIndex: Int = 0,
    /** 实际参与积分的窗口时长（秒），含顺延。 */
    val durationS: Double = 0.0,
    /** 按下到松手的时长（秒），即用户按住的时间。 */
    val holdS: Double = 0.0,
    val biasBound: Double = 0.0,
    val gravityMag: Double = 0.0,
    val accSigma: Double = 0.0,
    val dt: Double = 0.0,
) {
    val distanceCm: Double get() = distanceM * 100.0
    val rawDistanceCm: Double get() = rawDistanceM * 100.0
}

// ------------------------------------------------------------------ 入口

object SlideMeasure {

    /**
     * 一次测量。[downNs] / [upNs] 是按下与松手时刻的 `elapsedRealtimeNanos`，
     * 与 `Sample.tNs` 同基准，可直接二分对齐。
     *
     * [gScale] 保留只为与 Python 逐位对齐；**默认 1.0 = 不修正**，理由见文件头。
     */
    fun measure(
        acc: List<Pt3>,
        gyr: List<Pt3>,
        downNs: Long,
        upNs: Long,
        gScale: Double = 1.0,
    ): MeasureResult {
        val n = min(acc.size, gyr.size)
        if (n < MIN_SAMPLES) return MeasureResult(false, "样本太少（$n 条）")

        // ---- 1. 时间基准 ----
        val dts = ArrayList<Double>(n - 1)
        for (i in 0 until n - 1) dts.add((acc[i + 1].tNs - acc[i].tNs).toDouble())
        val dt = median(dts) / 1e9
        if (dt <= 0.0 || dt > 0.5) return MeasureResult(false, "时间戳异常（dt=${dt}s）")

        // ---- 2. 边界 ----
        val lo = locate(acc, downNs, n)
        val hi = locate(acc, upNs, n)
        if (lo < 20) return MeasureResult(false, "按下太早（校准段不足，lo=$lo）")
        if (hi <= lo + 5) return MeasureResult(false, "按钮窗口无效（lo=$lo hi=$hi）")
        if (hi >= n - 1) return MeasureResult(false, "松手太晚（尾部无数据）")

        // ---- 3. 自校准（DOWN 之前的静止段）----
        val calN = max(20, min((CAL_MAX_S / dt).toInt(), lo))
        val calLo = (lo - calN).coerceAtLeast(0)
        val cn = lo - calLo
        if (cn < 20) return MeasureResult(false, "静止段不足（$cn 条）")

        var gx = 0.0
        var gy = 0.0
        var gz = 0.0
        var ax = 0.0
        var ay = 0.0
        var az = 0.0
        for (i in calLo until lo) {
            gx += gyr[i].x; gy += gyr[i].y; gz += gyr[i].z
            ax += acc[i].x; ay += acc[i].y; az += acc[i].z
        }
        val inv = 1.0 / cn
        val gyroBias = V3(gx * inv, gy * inv, gz * inv)
        val aMean = V3(ax * inv, ay * inv, az * inv)
        val gravityRef = aMean.unit()
        val gravityMag = aMean.norm()
        if (gravityMag < 1.0) return MeasureResult(false, "重力读数异常（|g|=${gravityMag}）")

        // 噪声地板：逐轴取最大（不能三轴混算 —— z 轴含 9.8）
        val accSigma = maxOf(
            pstdev(calLo, lo) { acc[it].x },
            pstdev(calLo, lo) { acc[it].y },
            pstdev(calLo, lo) { acc[it].z },
        )
        val gyroSigma = maxOf(
            pstdev(calLo, lo) { gyr[it].x },
            pstdev(calLo, lo) { gyr[it].y },
            pstdev(calLo, lo) { gyr[it].z },
        )
        val profile = DeviceProfile(dt, gyroBias, gravityRef, gravityMag, accSigma, gyroSigma, calN)

        // ---- 4. 姿态：四元数互补滤波 ----
        val nq = n - 1
        val qs = Array(nq + 1) { Quat.IDENTITY }
        var q = qFromTwoVectors(profile.gravityRef, V3.UP)
        qs[0] = q
        for (i in 0 until nq) {
            val w = V3(gyr[i].x, gyr[i].y, gyr[i].z) - profile.gyroBias
            val qPred = q.mul(qDelta(w, dt)).normalized()
            val gPredBody = qPred.rotateInv(V3.UP)
            val e = gPredBody.cross(profile.gravityRef)
            q = if (e.norm() > REASON_OK_RAD) {
                qPred                                   // 振动防护：只外推
            } else {
                qPred.mul(qDelta(e * -ALPHA, dt)).normalized()
            }
            qs[i + 1] = q
        }

        // ---- 5. 去重力（用实测 |g|，不是 9.81）----
        val aLin = Array(n) { V3.ZERO }
        for (i in 0 until n) {
            val fw = qs[i].rotate(V3(acc[i].x, acc[i].y, acc[i].z))
            aLin[i] = V3(fw.x, fw.y, fw.z - profile.gravityMag)
        }

        // ---- 6~8. 距离 ----
        val calRangeLo = max(0, lo - (CAL_MAX_S / dt).toInt())
        val biasBound = meanVec(aLin, calRangeLo, lo).norm()

        val first = condZupt3d(aLin, lo, hi, dt, biasBound, accSigma, K_SIGMA, gScale)
        var end = hi
        if (abs(first.vEndFwd) > first.threshold) {
            end = stopEndpoint(aLin, lo, hi, dt, first.threshold, MAX_EXT_S, DWELL_S)
        }
        val final = if (end == hi) first
        else condZupt3d(aLin, lo, end, dt, biasBound, accSigma, K_SIGMA, gScale)

        val rawPath = integrateVec(aLin, lo, end, dt)
        val rawMag = rawPath.second[end - lo].norm()

        return MeasureResult(
            ok = true,
            distanceM = final.distanceM,
            rawDistanceM = rawMag,
            vEndFwd = final.vEndFwd,
            threshold = final.threshold,
            detrended = final.detrended,
            extended = end != hi,
            extendedMs = (end - hi) * dt * 1000.0,
            downIndex = lo,
            upIndex = hi,
            endIndex = end,
            durationS = (end - lo) * dt,
            holdS = (hi - lo) * dt,
            biasBound = biasBound,
            gravityMag = profile.gravityMag,
            accSigma = accSigma,
            dt = dt,
        )
    }
}

// ------------------------------------------------------------------ 内部类型

private class V3(val x: Double, val y: Double, val z: Double) {
    operator fun plus(o: V3) = V3(x + o.x, y + o.y, z + o.z)
    operator fun minus(o: V3) = V3(x - o.x, y - o.y, z - o.z)
    operator fun times(s: Double) = V3(x * s, y * s, z * s)
    fun dot(o: V3): Double = x * o.x + y * o.y + z * o.z
    fun cross(o: V3) = V3(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x)
    fun norm(): Double = sqrt(dot(this))
    fun unit(): V3 {
        val n = norm()
        return if (n > 1e-12) V3(x / n, y / n, z / n) else V3.UP
    }

    companion object {
        val ZERO = V3(0.0, 0.0, 0.0)
        val UP = V3(0.0, 0.0, 1.0)
    }
}

private class Quat(val w: Double, val x: Double, val y: Double, val z: Double) {
    fun mul(o: Quat) = Quat(
        w * o.w - x * o.x - y * o.y - z * o.z,
        w * o.x + x * o.w + y * o.z - z * o.y,
        w * o.y - x * o.z + y * o.w + z * o.x,
        w * o.z + x * o.y - y * o.x + z * o.w,
    )

    fun conj() = Quat(w, -x, -y, -z)

    fun normalized(): Quat {
        val n = sqrt(w * w + x * x + y * y + z * z)
        return if (n > 1e-12) Quat(w / n, x / n, y / n, z / n) else IDENTITY
    }

    /** body → world */
    fun rotate(v: V3): V3 {
        val r = mul(Quat(0.0, v.x, v.y, v.z)).mul(conj())
        return V3(r.x, r.y, r.z)
    }

    /** world → body */
    fun rotateInv(v: V3): V3 = conj().rotate(v)

    companion object {
        val IDENTITY = Quat(1.0, 0.0, 0.0, 0.0)
    }
}

private class DeviceProfile(
    val dt: Double,
    val gyroBias: V3,
    val gravityRef: V3,
    val gravityMag: Double,
    val accSigma: Double,
    val gyroSigma: Double,
    val nStatic: Int,
)

// ------------------------------------------------------------------ 算法原语

/** 角速度 [w] 持续 [dt] 的小旋转四元数。 */
private fun qDelta(w: V3, dt: Double) =
    Quat(1.0, w.x * dt / 2.0, w.y * dt / 2.0, w.z * dt / 2.0).normalized()

/** 求 q，使 R(q)·[aBody] == [bWorld]（两向量都当单位向量用）。 */
private fun qFromTwoVectors(aBody: V3, bWorld: V3): Quat {
    val a = aBody.unit()
    val b = bWorld.unit()
    val c = a.cross(b)
    val d = a.dot(b)
    if (d < -0.999999) {
        // 反向：绕任意垂直轴转 180°
        var axis = a.cross(V3(1.0, 0.0, 0.0))
        if (axis.norm() < 1e-6) axis = a.cross(V3(0.0, 1.0, 0.0))
        axis = axis.unit()
        return Quat(0.0, axis.x, axis.y, axis.z)
    }
    return Quat(1.0 + d, c.x, c.y, c.z).normalized()
}

/** 三维梯形二次积分，返回 (v, p)，下标 0..hi-lo。 */
private fun integrateVec(a: Array<V3>, lo: Int, hi: Int, dt: Double): Pair<Array<V3>, Array<V3>> {
    val m = hi - lo
    val v = Array(m + 1) { V3.ZERO }
    val p = Array(m + 1) { V3.ZERO }
    for (k in 1..m) {
        val i = lo + k
        v[k] = v[k - 1] + (a[i - 1] + a[i]) * 0.5 * dt
        p[k] = p[k - 1] + (v[k - 1] + v[k]) * 0.5 * dt
    }
    return v to p
}

private class Zupt(
    val distanceM: Double,
    val vEndFwd: Double,
    val threshold: Double,
    val detrended: Boolean,
)

/**
 * 三维净位移 + 条件矢量 ZUPT。
 *
 * 判据用「沿前进方向的终点速度」而不是三维模长 —— 模长里混着旋转/侧向的伪速度，
 * 会把阈值轻易顶穿、误触发去趋势。前进方向取 |v| 最大处的速度方向（那一刻最干净）。
 *
 * |v_end| ≤ thr → 手机已经停稳，零偏那点残余不值得动它 → 裸积分
 * |v_end| > thr → 显然还在动 → 三轴各自扣掉线性漂移（矢量 ZUPT）
 */
private fun condZupt3d(
    a: Array<V3>,
    lo: Int,
    hi: Int,
    dt: Double,
    biasBound: Double,
    sigma: Double,
    kSigma: Double,
    gScale: Double,
): Zupt {
    val m = hi - lo
    val (v, p) = integrateVec(a, lo, hi, dt)
    val t = m * dt

    var kpk = 1
    var bestNorm = -1.0
    for (k in 1..m) {
        val nn = v[k].norm()
        if (nn > bestNorm) { bestNorm = nn; kpk = k }
    }
    val u = v[kpk].unit()
    val vEndFwd = v[m].dot(u)
    val thr = biasBound * t + kSigma * sigma * sqrt(max(t * dt, 0.0))

    if (abs(vEndFwd) <= thr) {
        return Zupt(p[m].norm() / gScale, vEndFwd, thr, false)
    }
    val dr = if (t > 1e-9) V3(v[m].x / t, v[m].y / t, v[m].z / t) else V3.ZERO
    val a2 = Array(m + 1) { i -> a[lo + i] - dr }
    val p2 = integrateVec(a2, 0, m, dt).second
    return Zupt(p2[m].norm() / gScale, vEndFwd, thr, true)
}

/**
 * 把终点从「松手点 hi」顺延到「手机真正停稳那一刻」。
 *
 * 只在松手时手机还在动时调用。不能用「速度过零」找停稳点 —— 从 hi 起积分速度天生为 0；
 * 这里的速度一律从 lo（按下点）积分，才是绝对速度。
 *
 *   1. 找第一个「|v_fwd| ≤ thr 且持续 dwell」的样本 → 停稳区起点 kStop
 *   2. 从 kStop 往后扩到 |v_fwd| 重新超阈为止 → 停稳区终点 kEnd
 *   3. 终点取 [hi, kEnd] 内**前向位移最大**的那一帧：手机停住后可能原地不动（平台），
 *      也可能被手拉回去（折返）—— 两种情况的"位移见顶"都指向同一个正确终点。
 *
 * 找不到合法停稳区就原样返回 hi。
 */
private fun stopEndpoint(
    a: Array<V3>,
    lo: Int,
    hi: Int,
    dt: Double,
    thr: Double,
    maxExtS: Double,
    dwellS: Double,
): Int {
    val lim = min(a.size - 1, hi + (maxExtS / dt).toInt())
    if (lim - hi < 5) return hi

    val (v, p) = integrateVec(a, lo, lim, dt)
    val kHi = hi - lo
    val u = (p[kHi] - p[0]).unit()
    val fv = DoubleArray(v.size) { v[it].dot(u) }
    val fp = DoubleArray(p.size) { p[it].dot(u) }
    val dwell = max(2, (dwellS / dt).toInt())

    var kStop = -1
    for (k in kHi + 1 until fv.size) {
        var settled = true
        for (j in k until min(fv.size, k + dwell)) {
            if (abs(fv[j]) > thr) { settled = false; break }
        }
        if (settled) { kStop = k; break }
    }
    if (kStop < 0) return hi

    var kEnd = kStop
    for (k in kStop until fv.size) {
        if (abs(fv[k]) <= thr) kEnd = k else break
    }
    if (kEnd - kStop < dwell) return hi

    var best = kHi
    for (k in kHi..kEnd) if (fp[k] > fp[best]) best = k
    return if (best > kHi) lo + best else hi
}

// ------------------------------------------------------------------ 小工具

/** 二分找第一个 `tNs >= target` 的下标（[limit] 是可用长度）。 */
private fun locate(xs: List<Pt3>, target: Long, limit: Int): Int {
    var lo = 0
    var hi = limit - 1
    while (lo < hi) {
        val mid = (lo + hi) / 2
        if (xs[mid].tNs < target) lo = mid + 1 else hi = mid
    }
    return lo
}

/** 与 Python `statistics.median` 一致：偶数个取中间两个的平均。 */
private fun median(xs: List<Double>): Double {
    if (xs.isEmpty()) return 0.0
    val s = xs.sorted()
    val n = s.size
    return if (n % 2 == 1) s[n / 2] else (s[n / 2 - 1] + s[n / 2]) / 2.0
}

/** 总体标准差（除以 N，与 Python `statistics.pstdev` 一致）。 */
private inline fun pstdev(from: Int, to: Int, pick: (Int) -> Double): Double {
    val n = to - from
    if (n <= 0) return 0.0
    var sum = 0.0
    for (i in from until to) sum += pick(i)
    val m = sum / n
    var acc = 0.0
    for (i in from until to) {
        val d = pick(i) - m
        acc += d * d
    }
    return sqrt(acc / n)
}

private fun meanVec(xs: Array<V3>, from: Int, to: Int): V3 {
    val n = to - from
    if (n <= 0) return V3.ZERO
    var x = 0.0
    var y = 0.0
    var z = 0.0
    for (i in from until to) {
        x += xs[i].x; y += xs[i].y; z += xs[i].z
    }
    return V3(x / n, y / n, z / n)
}
