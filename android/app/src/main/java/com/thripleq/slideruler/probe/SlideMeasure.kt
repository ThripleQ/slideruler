package com.thripleq.slideruler.probe

import java.util.Locale
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt

/**
 * 测距引擎 v2：先把「静止」和「运动」分开，再在两个静止锚点之间积分。
 *
 * ## v1 错在哪（用户报「按着不动也能测出很长一段距离」）
 *
 * v1 用「条件 ZUPT」：`|v_end| ≤ 阈值` 就走裸积分。这个判据**自相矛盾** —— 它拿零偏估计
 * `b` 去判断"终点速度是不是零偏造成的"，而零偏正是误差来源。于是任何静止窗口都必然通过
 * 检验（`|v_end| = b·T` 恒等于阈值的第一项），必然输出 `½·b·T²` 的假位移。
 *
 * 真机实测：`20261011_100631`，手机 3.64 s 全程静止（窗口内静止样本 100%），
 * `|b|=0.035 m/s²` → 裸积分 **18.20 cm**。
 *
 * ## v2 的三条原则
 *
 * 1. **静止检测只用原始数据**（加速度局部方差 + 陀螺模长），不依赖姿态、不依赖零偏。
 *    这样它不会和它要保护的东西互相污染。
 * 2. **窗口 = 真实运动段**，由前后两段验证过的静止夹出来 —— 不是用户按住的时长。
 *    静止段里的真实位移恒为 0，却贡献完整的零偏误差（∝T²），所以多算一秒纯亏。
 *    用户「按住→停顿→推→停稳→松手」里的停顿和停稳，本该被切掉。
 * 3. **两端都验证过静止 ⇒ 无条件去趋势**，而且把两侧静止段里的**所有**样本一起当锚点
 *    （双锚点 ZUPT），而不是只钉两个端点。噪声按 √N 平均掉。
 *
 * ## 与 v1 的其它差别
 *
 * · 姿态：加速度修正的参考从「校准时刻的体重力方向」（体坐标里固定）换成**当前实测**
 *   加速度方向（标准 Mahony），且只在静止段开修正。v1 的写法等于让滤波器主动抵抗手机的
 *   真实旋转，持续角速率 ω 下稳态姿态误差 ≈ ω·dt/α（ω=5°/s 就有 1°、重力泄漏 0.17 m/s²）。
 * · 校准：从「DOWN 之前 2 s」换成**运动之前那段验证过的静止**（就在按钮窗内，用户的停顿）。
 *   实测 20261011_100339 就是被这条坑的：用户按下前正在挪手机，写死的 2 s 里 |g|=9.43
 *   （应该 9.91）、σ_a=2.01 → 报告 168 cm。
 * · 拿不到两端静止就**拒答**，并说清为什么 —— 不给一个看起来体面的假数字。
 * · 「平台宽度」= 把窗口端点在各侧静止段内挪动，距离的变化范围。这是**自带的误差棒**：
 *   平台窄 = 判断稳；平台宽 = 换个端点就变，这数不可信。
 *
 * 每个数字的推导见 `sim/measure2.py`（Python 参考实现，两者逐位对齐）与
 * `sim/exp_v2.py`（有精确真值的仿真对照）。
 */

// ------------------------------------------------------------------ 参数

private const val STATIC_WIN_S = 0.20      // 静止检测滑动窗（边界要锐）
private const val SLOW_WIN_S = 0.80        // 慢动作检测窗（白噪声 σ 不随窗长变，慢信号会涨）
private const val STATIC_FACTOR = 4.0      // 静止阈值 = 本机噪声地板 × 该系数
private const val STATIC_FLOOR = 0.010     // 噪声地板硬下限（m/s²）
private const val GYRO_FLOOR = 0.010       // 陀螺静止阈值硬下限（rad/s）
private const val STATIC_MIN_S = 0.20      // 一段静止至少这么久才敢当锚点（约束整段，见 measure()）
private const val CAL_MIN_S = 0.25         // 校准段最短
private const val CAL_MAX_S = 2.0          // 校准段最长
private const val CAL_STATIC_FRAC = 0.95   // 校准段里静止样本的最低占比（低于此拒绝测量）
private const val MARGIN_S = 0.15          // 运动段两侧各留的静止余量
private const val TAIL_S = 2.5             // UP 之后最多再看这么久
private const val PRE_S = 2.0              // DOWN 之前最多往回看这么久
private const val GAP_S = 0.15             // 两段运动之间隔这么久以内算同一段
private const val ALPHA = 0.05             // 互补滤波系数
private const val GATE_RAD = 0.30          // 加速度修正门限（预测/实测重力方向夹角），仅 gMag 未知时用
private const val GATE_MAG_FRAC = 0.25     // |acc| 偏离 |g| 超过这个比例就不修正（见 solveAttitude）
private const val PLATEAU_GRID = 5         // 平台扫描点数
// 「修正量相对答案」的上限：|裸积分 − 去趋势| 超过 max(答案, 该下限) 就报"模型敏感"。
// 0.02 m 这个下限是为了避免小距离时被绝对噪声误报。
private const val MODEL_GAP_FLOOR_M = 0.02
// 运动段速度峰值低于此值时，「锚点速度 / 运动峰值」是两个噪声相除，判据无意义
// （真机 122256 的运动峰值只有 0.003 m/s，几乎没动，比值却到 1.03）。
private const val ANCHOR_MIN_PEAK_MPS = 0.05
private const val MIN_SAMPLES = 40

/** 一条三轴样本。用 Double 而不是 Float，是为了让单测能与 Python 逐位对齐。 */
class Pt3(val tNs: Long, val x: Double, val y: Double, val z: Double)

fun Sample.toPt3(): Pt3 = Pt3(tNs, x.toDouble(), y.toDouble(), z.toDouble())

/**
 * 测量结果。[ok] 为 false 时 [failReason] 给出原因。
 *
 * 除了距离，刻意把「引擎看到的证据」一并带出来 —— 用户需要知道这一次可不可信：
 * [confidence] 是自带的误差棒，[noMotion]/[slowMotion]/[endHitBound] 是三种"没测好"的具体原因。
 *
 * `data class` 是为了 `copy()` —— 模型敏感时要把同一个结果翻成"拒答"，
 * 又不想把三十来个字段再抄一遍（抄一遍就会与成功路径漂移）。
 */
data class MeasureResult(
    val ok: Boolean,
    val failReason: String? = null,
    /** **测量头**的净位移（米）。headR 全 0 时等于 [imuDistM]。 */
    val distanceM: Double = 0.0,
    /** IMU 芯片自身的净位移（米）。与 [distanceM] 的差就是杠杆臂项 [leverM] 的作用。 */
    val imuDistM: Double = 0.0,
    /** 杠杆臂修正量 |[R(t1)−R(t0)]·r|（米）。比读数还大 ⇒ 姿态变化主导了这次测量。 */
    val leverM: Double = 0.0,
    /** 同一窗口**不去趋势**（裸积分）得到的位移（米）。与 [distanceM] 的差 = 修正量。 */
    val bareM: Double = 0.0,
    /** |裸积分 − 去趋势|（米）：残差模型的分歧有多大。 */
    val modelGapM: Double = 0.0,
    /** 修正量超过答案本身 —— 这个数由**模型假设**决定、不是由数据决定。 */
    val modelSensitive: Boolean = false,
    /** 平台下/上界（米）：把窗口端点在各侧静止段内挪动时距离的变化范围。 */
    val plateauLoM: Double = 0.0,
    val plateauHiM: Double = 0.0,
    /** 窗口内确实没有运动（距离按 0 处理）。 */
    val noMotion: Boolean = false,
    /** 有运动但太慢，检测器分辨不出（此时不给距离）。 */
    val slowMotion: Boolean = false,
    /** 真实运动时长（秒）—— 已被裁掉停顿与停稳。 */
    val motionS: Double = 0.0,
    /** 实际积分窗口时长（秒）。 */
    val windowS: Double = 0.0,
    /** 用户按住时长（秒）。 */
    val holdS: Double = 0.0,
    /** 运动段结束后没在 [TAIL_S] 内找到静止段 —— 松手时手机还在动，或尾部锚点贴在数据末尾。 */
    val endHitBound: Boolean = false,
    /** 松手之后到底留了多少静止当终点锚点（秒）。< 0.15 s 就是"没得用"，这个数不可靠。 */
    val tailStaticS: Double = 0.0,
    /** 修正前速度在两侧静止段里的均值模长（米/秒）—— 也就是"要去掉多少速度"。 */
    val anchorDrift: Double = 0.0,
    val motionIndex: Int = 0,
    /**
     * 这次测量合并了几段被静止检测切开的运动。
     * 1 = 一次连续滑动；>1 = 中间有过停顿或匀速段被误判成静止，已合并。
     */
    val mergedRuns: Int = 1,
    val windowStart: Int = 0,
    val windowEnd: Int = 0,
    val anchorStart: Int = 0,
    val anchorEnd: Int = 0,
    val calS: Double = 0.0,
    val gravityMag: Double = 0.0,
    val accSigma: Double = 0.0,
    val gyroBiasNorm: Double = 0.0,
    val dt: Double = 0.0,
    val slowPeak: Double = 0.0,
    val slowThr: Double = 0.0,
    /** 陀螺零偏**不稳定性**（rad/s，τ 尺度）：静止段内分段均值之差 /√2。抹不掉的那部分。 */
    val gyroInstability: Double = 0.0,
    /** 上面的估计所依据的时间尺度（秒）。 */
    val gyroInstTauS: Double = 0.0,
    /**
     * 姿态漂移造成的位移误差**上界**（米）：(1/6)·g·σ_bg·√(T/τ)·T³。
     * **∝ T³** —— 滑动时间减半，误差降到 1/8，比任何算法微调都管用。
     * 保守上界（假设游走全朝一个方向），只作操作引导，不据此拒答。
     */
    val attBudgetM: Double = 0.0,
    /** 窗口内**运动段**上裸速度的峰值（米/秒）。只用来判断终点锚点可信度。 */
    val vPeakMps: Double = 0.0,
    /**
     * 终点锚点可疑：标为「静止」的那段里，积分速度与运动段速度峰值同量级。
     *
     * 「松手」不等于「停住」，而**匀速移动与静止在加速度计上完全不可分辨**（测的是比力，
     * 匀速时比力就是重力）—— 静止检测必然把匀速段判成静止。双锚点于是照假设把它当零，
     * 把**真实位移**当漂移扣掉：真机 20261011_143734（真值 15 cm）终点锚点段速度
     * 0.316 m/s、运动段峰值 0.314 m/s，扣掉 15.42 cm，16.19 → 1.59 cm。
     * v3 那条教训（匀速段把一次滑动切碎）当时只用来改运动段判据，没意识到它同样毁掉锚点。
     *
     * **只诊断、不拒答**：143725 同样可疑（0.90）却给出 +5% 的正确读数（两个错误互相抵消）。
     * 这条判据没有干净的分界线，拿它拒答等于拿真值调参 —— 本项目不做。
     */
    val anchorSuspect: Boolean = false,
) {
    val distanceCm: Double get() = distanceM * 100.0
    val plateauLoCm: Double get() = plateauLoM * 100.0
    val plateauHiCm: Double get() = plateauHiM * 100.0
    /** 误差棒 / 读数。0.02 以内算稳，0.1 以上基本不能用。 */
    val confidence: Double
        get() = if (distanceM > 1e-9) (plateauHiM - plateauLoM) / distanceM else 0.0
}

// ------------------------------------------------------------------ 入口

object SlideMeasure {

    /** 平台扫描 / 慢动作守卫用的参数，暴露给单测。 */
    const val TAIL_SECONDS = TAIL_S

    /** 终点锚点至少要这么多静止才算数（= 窗口余量 MARGIN_S）。 */
    const val TAIL_MIN_SECONDS = MARGIN_S

    /**
     * 一次测量。[downNs] / [upNs] 是按下与松手时刻的 `elapsedRealtimeNanos`，
     * 与 `Sample.tNs` 同基准，可直接二分对齐。
     */
    fun measure(
        acc: List<Pt3>, gyr: List<Pt3>, downNs: Long, upNs: Long,
        /**
         * IMU → 测量头的机体系偏移（米）。默认全 0 = 把 IMU 自己当测量头（v3 行为）。
         * 非 0 时输出的是**测量头**的净位移：Δp_H = Δp_I + [R(t1) − R(t0)]·r
         * —— 刚体关系在位置层面就是精确的，不需要对角速度求差分。
         * ⚠ 这个向量必须标定（见 docs/v4-head-model.md §3），猜出来的值会让修正反向。
         */
        headRx: Double = 0.0, headRy: Double = 0.0, headRz: Double = 0.0,
    ): MeasureResult {
        val n = min(acc.size, gyr.size)
        if (n < MIN_SAMPLES) return MeasureResult(false, "样本太少（$n 条）")

        // ---- 1. 时间基准 ----
        var dt = 0.0
        run {
            val dts = DoubleArray(n - 1)
            for (i in 0 until n - 1) dts[i] = (acc[i + 1].tNs - acc[i].tNs).toDouble()
            dt = median(dts) / 1e9
        }
        if (dt <= 0.0 || dt > 0.5) return MeasureResult(false, "时间戳异常（dt=${dt}s）")

        val ax = DoubleArray(n) { acc[it].x }
        val ay = DoubleArray(n) { acc[it].y }
        val az = DoubleArray(n) { acc[it].z }
        val gx = DoubleArray(n) { gyr[it].x }
        val gy = DoubleArray(n) { gyr[it].y }
        val gz = DoubleArray(n) { gyr[it].z }
        val dn = locate(acc, downNs, n)
        val up = locate(acc, upNs, n)
        // 按钮时段必须有正长度。DOWN==UP 时所有运动组与按钮窗的重叠都是 0，
        // 下面「取重叠最大的一组」就退化成「取时间最早的一组」—— 会去测**按住之前**
        // 那一小段（用户挪手机）而不是这次滑动。实测 100444 上会报 0.06 cm 而不是 15.4 cm。
        if (up <= dn) return MeasureResult(false, "按下与松手时刻异常（松手不晚于按下）")

        // ---- 2. 静止检测（只用原始数据）----
        //      陀螺零偏先用「最安静的 0.5 s」粗估，只为了把 ω 的直流分量去掉，精度要求很低。
        val sdFast = rollingMaxStd(ax, ay, az, dt, STATIC_WIN_S)
        var iQuiet = 0
        for (i in 1 until n) if (sdFast[i] < sdFast[iQuiet]) iQuiet = i
        val halfPre = max(10, min(n / 4, (0.5 / dt).toInt()) / 2)
        val qa = max(0, iQuiet - halfPre)
        val qb = min(n, iQuiet + halfPre + 1)
        var gbx = 0.0
        var gby = 0.0
        var gbz = 0.0
        for (i in qa until qb) { gbx += gx[i]; gby += gy[i]; gbz += gz[i] }
        val invQ = 1.0 / (qb - qa)
        gbx *= invQ; gby *= invQ; gbz *= invQ

        val floorA = percentile(sdFast, 0.05)
        val thrA = max(STATIC_FACTOR * floorA, STATIC_FLOOR)
        val gn = DoubleArray(n) {
            val ex = gx[it] - gbx; val ey = gy[it] - gby; val ez = gz[it] - gbz
            sqrt(ex * ex + ey * ey + ez * ez)
        }
        val gnRoll = rollingMax(gn, dt, STATIC_WIN_S)
        val floorG = percentile(gnRoll, 0.05)
        val thrG = max(STATIC_FACTOR * floorG, GYRO_FLOOR)
        val flags = BooleanArray(n) { sdFast[it] <= thrA && gnRoll[it] <= thrG }

        val minRun = max(3, (STATIC_MIN_S / dt).toInt())

        // ---- 3. 搜索区：DOWN 前 2 s → UP 后 TAIL_S ----
        val regLo = max(0, dn - (PRE_S / dt).toInt())
        val regHi = min(n - 1, up + (TAIL_S / dt).toInt())
        val gap = max(2, (GAP_S / dt).toInt())

        val mot = ArrayList<Int>()
        for (i in regLo..regHi) if (!flags[i]) mot.add(i)
        if (mot.isEmpty()) {
            // 0.2 s 的局部标准差看不出一场很慢的滑动（15 cm 滑 4 s，峰值加速度只有
            // 0.059 m/s²，0.2 s 窗内的变化量 ≈0.006 —— 埋在噪声里）。0.8 s 尺度上看得见。
            // 这里不猜距离，只把"其实移动了"说出来，让用户滑快一点。
            val sdSlow = rollingMaxStd(ax, ay, az, dt, SLOW_WIN_S)
            var peak = 0.0
            for (i in regLo..regHi) if (sdSlow[i] > peak) peak = sdSlow[i]
            val slowThr = max(STATIC_FACTOR * percentile(sdSlow, 0.05), STATIC_FLOOR)
            if (peak > slowThr) {
                return MeasureResult(false,
                    "移动太慢，检测器分辨不出 —— 请快速滑完（1~2 秒）",
                    slowMotion = true, slowPeak = peak, slowThr = slowThr, dt = dt)
            }
            return MeasureResult(true, distanceM = 0.0, noMotion = true, dt = dt,
                holdS = (up - dn) * dt, slowPeak = peak, slowThr = slowThr)
        }

        // 按间隙分组
        val groups = ArrayList<IntArray>()
        var cur = intArrayOf(mot[0], mot[0])
        for (k in 1 until mot.size) {
            val i = mot[k]
            if (i - cur[1] <= gap) cur[1] = i else { groups.add(cur); cur = intArrayOf(i, i) }
        }
        groups.add(cur)

        // ---- 运动区间的选法：按钮按住期间发生过的运动，取并集 ----
        // 不能靠静止检测去切运动段。**匀速运动在加速度计上和静止无法区分**
        // （加速度计测比力；匀速时比力就是重力，和静止完全一样），
        // 所以一次平稳滑动的中间段必然被判成"静止"，完整滑动被切成好几段。
        // 真机证据 122416：一次按住里被切成 5 段（0.20/0.29/0.14/0.26/1.19 s），
        // 只取其中一段 → 报 0.59 cm（看起来体面但错）。
        //
        // 旧的「取与按钮窗重叠最多的那一组」有两个毛病：
        //   ① 所有组都不与按钮窗重叠时退化成「取时间最早的一组」—— 那是按下前挪手机，
        //      不是这次滑动（`up <= dn` 那个守卫挡不住：它只保证按钮窗长度为正）。
        //   ② **惩罚延续到松手之后的那段**：正常滑动常常滑过松手时刻，它在窗内的重叠
        //      反而小于一段完全落在窗内的短促抖动。实测 122416：真正那次 1.19 s 的滑动
        //      重叠只有 0.24 s，输给了 0.29 s 的抖动段（重叠 0.29 s）—— 选中的是错的。
        var m0 = -1
        var m1 = -1
        var nSel = 0
        for (g in groups) {
            if (max(0, min(g[1], up) - max(g[0], dn)) <= 0) continue
            nSel++
            if (m0 < 0 || g[0] < m0) m0 = g[0]
            if (g[1] > m1) m1 = g[1]
        }
        if (nSel == 0) {
            return MeasureResult(false,
                "按下期间没有检测到运动（${groups.size} 段运动都在按下之前或松手之后）",
                dt = dt, holdS = (up - dn) * dt)
        }

        // ---- 4. 窗口 = 运动段 ± 余量，两端必须落在静止段里 ----
        val marg = max(1, (MARGIN_S / dt).toInt())
        var w0 = m0 - marg
        var w1 = m1 + marg
        while (w0 > 0 && !flags[w0]) w0--
        while (w1 < n - 1 && !flags[w1]) w1++
        w0 = max(0, w0)
        w1 = min(n - 1, w1)
        // 「终点锚点没得用」有两种：① 运动一直延续到搜索区末尾（手机没停稳）；
        // ② 算法要把窗口往右放 MARGIN_S 的余量，但数据在 w1=n-1 处就到头了 ——
        //    这时终点锚点会被裁成只剩几个样本甚至 1 个，白噪声平均不掉。
        //    实测 100555 尾部静止只有 0.047 s，终点锚点退化成 1 个样本，
        //    它给出的 1.01 cm 完全靠不住，却只由平台宽度（60%）间接暴露。
        val endHitBound = m1 >= regHi - 1 || w1 >= n - 1
        if (!flags[w0] || !flags[w1]) {
            // 两端要分开说。写死一句「松手后没停稳」会让用户朝着错误方向修 ——
            // 实测 100601 / 100608 是**按下之前**手机一直在动（搜索区里没有一个静止样本），
            // 用户再怎么「松手后保持不动」也没用，得先放稳再按。
            val parts = buildList {
                if (!flags[w0]) add("按下之前手机一直在动")
                if (!flags[w1]) add("松手后手机没停稳")
            }
            return MeasureResult(false,
                "找不到运动两端的静止段（${parts.joinToString("，")}）",
                endHitBound = endHitBound, motionIndex = m0, dt = dt, holdS = (up - dn) * dt,
                motionS = (m1 - m0) * dt, windowStart = w0, windowEnd = w1)
        }
        if (w1 - w0 < max(5, (0.05 / dt).toInt())) {
            return MeasureResult(false, "有效窗口太短", dt = dt)
        }

        // 窗口两端的整段静止（ZUPT 锚点；在 zupt 内会裁到窗口，所以窗口不会因此变长）
        val run0 = runContaining(flags, w0, n)
        val run1 = runContaining(flags, w1, n)
        if (run0[0] == run1[0]) {
            return MeasureResult(true, distanceM = 0.0, noMotion = true, dt = dt,
                holdS = (up - dn) * dt)
        }
        // STATIC_MIN_S 约束的是**整段静止**，不是裁进窗口之后剩下的那点余量。
        // 后者恒等于 MARGIN_S（实测 0.143 s / 15 个样本，见 sim/exp_anchor.py）——
        // 拿它去卡 0.20 s 会把全部测量干掉，那是约束错了对象。
        //
        // 想让有效锚点变长只有一个办法：把窗口往静止段里外扩。试过了，不行 ——
        // 三条有真值的里两条变差（100459 +5.02%→+9.37%、100529 +3.72%→+6.04%），
        // 偏差随窗口长度单调增长。静止段延长窗口**不是**白送样本（零偏会漂，
        // 线性模型拟合不住），∝T² 那一派赢了。所以锚点就这么长，认了。
        //
        // 那 STATIC_MIN_S 剩下的用处是挡「假静止」：0.2 s 的滚动标准差根本判不出
        // 只有 0.038 s 的静止，而 100555 右侧静止段只有 4 个样本，它给出的
        // 1.01 cm（平台 0.90~1.51）纯属噪声。整段静止够长，才说明这里真静过。
        if (min(run0[1] - run0[0], run1[1] - run1[0]) < minRun) {
            val s0 = (run0[1] - run0[0]) * dt
            val s1 = (run1[1] - run1[0]) * dt
            val which = if (s0 <= s1) "按下之前" else "松手之后"
            return MeasureResult(false,
                "${which}的静止段只有 ${"%.2f".format(min(s0, s1))} s" +
                    "（要 ≥${"%.2f".format(STATIC_MIN_S)} s 才敢当锚点）—— 停稳后再松手",
                dt = dt, holdS = (up - dn) * dt, motionS = (m1 - m0) * dt,
                windowStart = w0, windowEnd = w1)
        }

        // ---- 5. 校准段 = 运动之前那段静止（用户「按住后停顿」的那一段）----
        val calLo = max(run0[0], m0 - (CAL_MAX_S / dt).toInt())
        val calHi = m0
        if ((calHi - calLo) * dt < CAL_MIN_S) {
            return MeasureResult(false,
                "按下后先停一下再推（运动之前的静止段只有 ${"%.2f".format((calHi - calLo) * dt)} s）",
                dt = dt)
        }
        var cgx = 0.0; var cgy = 0.0; var cgz = 0.0
        var cax = 0.0; var cay = 0.0; var caz = 0.0
        for (i in calLo until calHi) {
            cgx += gx[i]; cgy += gy[i]; cgz += gz[i]
            cax += ax[i]; cay += ay[i]; caz += az[i]
        }
        val invC = 1.0 / (calHi - calLo)
        val bias = V3(cgx * invC, cgy * invC, cgz * invC)
        val aMean = V3(cax * invC, cay * invC, caz * invC)
        val gMag = aMean.norm()
        val gRef = aMean.unit()
        val sigma = maxOf(
            pstdev(calLo, calHi) { ax[it] },
            pstdev(calLo, calHi) { ay[it] },
            pstdev(calLo, calHi) { az[it] },
        )
        if (gMag < 0.5 * 9.80665 || gMag > 1.5 * 9.80665) {
            return MeasureResult(false, "校准段重力读数异常 |g|=${"%.3f".format(gMag)}", dt = dt)
        }
        // 校准段必须**几乎全静**。混着运动样本的均值既不是静止零偏也不是运动加速度，
        // 拿它当重力方向会让整条链路的重力泄漏掉一大块。
        // 实测 100339：|g| 读到 9.43（本机实测应 9.91）、σ=2.01，v1 据此报了 168.62 cm。
        // 这一条 Python 有、Kotlin 原本没有 —— 是移植时漏掉的，夹具没覆盖所以单测没发现。
        var staticCnt = 0
        for (i in calLo until calHi) if (flags[i]) staticCnt++
        val calStaticFrac = staticCnt.toDouble() / (calHi - calLo)
        if (calStaticFrac < CAL_STATIC_FRAC) {
            return MeasureResult(false,
                "校准段不够静止（静止占比 ${"%.0f".format(calStaticFrac * 100)}%）",
                dt = dt, calS = (calHi - calLo) * dt, gravityMag = gMag)
        }

        // ---- 6. 姿态（加速度修正只在静止段开，运动中纯陀螺外推）----
        // 按钮窗内一律不做加速度修正：那段时间手机在动，acc = 重力 + 运动加速度，
        // 它的方向**不是**重力方向。只影响姿态，不影响锚点搜索 —— 锚点仍用原始静止判据找。
        //
        // ⚠ 这一条移植时漏过（Python 有、Kotlin 没有），后果是 v4 的关键改动在设备端失效：
        //   flags 与 gate_flags 平时看不出差别，因为运动段的 |acc| 偏离 |g| 很大、
        //   门控本来就会拒；但**匀速段**的 |‖acc‖−g| ≈ 0 **会被门控放行** ——
        //   而匀速段恰恰是静止检测必然误判成静止的地方（匀速与静止在加速度计上不可分辨）。
        //   于是那些样本被拿去"修正"姿态，把运动加速度写进了重力方向。
        //   真机 122416（一次按住切 5 段）两侧因此差 9.4 cm，而单测当时没覆盖这种夹具。
        val gateFlags = BooleanArray(n) { flags[it] && (it < dn || it > up) }
        val qs = solveAttitude(ax, ay, az, gx, gy, gz, bias, gRef, gateFlags, dt,
            i0 = calLo, gMag = gMag)

        // ---- 7. 去重力（用实测 |g|，不是 9.81）----
        val lx = DoubleArray(n); val ly = DoubleArray(n); val lz = DoubleArray(n)
        for (i in 0 until n) {
            val f = qs[i].rotate(V3(ax[i], ay[i], az[i]))
            lx[i] = f.x; ly[i] = f.y; lz[i] = f.z - gMag
        }

        // ---- 8. 双锚点 ZUPT ----
        val zupt = zuptTwoAnchor(lx, ly, lz, w0, w1, run0, run1, dt)

        // ---- 8a. 测量头几何（刚体杠杆臂）----
        // 加速度计装在机身某处，用户想测的是机身另一个点，两者**不是同一个点**：
        // 手机一转，即使测量头钉着不动，IMU 也会沿弧线走一段。
        //     Δp_H = Δp_I + [R(t1) − R(t0)]·r
        // 在**位置层面**做修正，就不需要角加速度 —— 若改在加速度层面补偿
        // （a_H = a_I + R(α×r + ω×(ω×r))），得先对带噪的陀螺求差分，噪声会被放大。
        // 量级：r=9 cm、转 30° 时 |ΔR·r| = 2r·sin(15°) ≈ 4.7 cm。仿真里 v3 就报了这个数（真值 0）。
        val dLever = qs[w1].rotate(V3(headRx, headRy, headRz)) -
            qs[w0].rotate(V3(headRx, headRy, headRz))
        val pHead = V3(zupt.endX, zupt.endY, zupt.endZ) + dLever
        val distHead = pHead.norm()

        // ---- 8b. 残差模型分歧 ----
        // 「裸积分」与「去趋势」是残差加速度的两种极端模型：前者当它是 0，后者当它是 A·t+B。
        // 两者相差超过答案本身 ⇒ 答案是两个大数之差，由**模型假设**决定而不是由数据决定，
        // 这时窄平台毫无意义（平台只反映端点抖动，反映不了模型分歧）。
        // 真机实测的分界很干净：有三条有真值（15 cm）的是 3% / 16% / 46%，都没超；
        // 四条没有真值的是 150% / 260% / 590% / 4587%，全超。见 sim/compare_v1_v2.py。
        // 裸积分也要加同一个杠杆臂项，否则两个模型在不同对象上比，gap 没有意义。
        val bv = integrateEndVec(lx, ly, lz, w0, w1, dt)
        val bare = (V3(bv[0], bv[1], bv[2]) + dLever).norm()
        val modelGap = abs(bare - distHead)
        val modelSensitive = modelGap > max(distHead, MODEL_GAP_FLOOR_M)

        // ---- 9. 平台宽度：窗口端点在各自静止段内挪动，答案变化多少 ----
        var plo = Double.MAX_VALUE
        var phi = -Double.MAX_VALUE
        for (k in 0 until PLATEAU_GRID) {
            val f = k.toDouble() / (PLATEAU_GRID - 1)
            val a0 = Math.round(run0[0] + f * (m0 - run0[0])).toInt()
            for (j in 0 until PLATEAU_GRID) {
                val g = j.toDouble() / (PLATEAU_GRID - 1)
                val b1 = Math.round(m1 + g * (run1[1] - 1 - m1)).toInt()
                if (b1 - a0 < 5) continue
                val z = zuptTwoAnchor(lx, ly, lz, a0, b1,
                    runContaining(flags, a0, n), runContaining(flags, b1, n), dt)
                // 端点在动 ⇒ 杠杆臂项跟着端点姿态走，必须一起扫，
                // 否则平台只反映 IMU 那部分，掩盖了姿态误差经由杠杆臂的放大。
                val dl = qs[b1].rotate(V3(headRx, headRy, headRz)) -
                    qs[a0].rotate(V3(headRx, headRy, headRz))
                val d = (V3(z.endX, z.endY, z.endZ) + dl).norm()
                if (d < plo) plo = d
                if (d > phi) phi = d
            }
        }
        if (plo > phi) { plo = distHead; phi = distHead }

        // ---- 10. 姿态漂移的误差上界 ∝ T³ ----
        // "对操作要求高"的物理根源：陀螺零偏的游走抹不掉，它在运动段里变成
        // 姿态漂移 → 重力泄漏 → 位移误差 ∝ T³。滑快一倍，误差降到 1/8。
        // 只作操作引导，不据此拒答 —— 它是保守上界（假设游走全朝一个方向）。
        val inst = gyroBiasInstability(gx, gy, gz, run0, dt)
        val attBudget = attitudeDriftBudget(gMag, inst[0], inst[1], (w1 - w0) * dt)

        // ---- 11. 终点锚点可信度（只诊断，不拒答；理由见 MeasureResult.anchorSuspect）----
        val vPeak = motionVelPeak(lx, ly, lz, w0, w1, m0, m1, dt)
        val anchorSuspect = vPeak > ANCHOR_MIN_PEAK_MPS && zupt.anchorDrift > 0.5 * vPeak

        val out = MeasureResult(
            ok = true,
            distanceM = distHead,
            imuDistM = zupt.distanceM,
            leverM = dLever.norm(),
            bareM = bare,
            modelGapM = modelGap,
            modelSensitive = modelSensitive,
            plateauLoM = plo,
            plateauHiM = phi,
            motionS = (m1 - m0) * dt,
            windowS = (w1 - w0) * dt,
            holdS = (up - dn) * dt,
            endHitBound = endHitBound,
            tailStaticS = (run1[1] - m1) * dt,
            anchorDrift = zupt.anchorDrift,
            motionIndex = m0,
            mergedRuns = nSel,
            windowStart = w0,
            windowEnd = w1,
            anchorStart = run0[0],
            anchorEnd = run1[1] - 1,
            calS = (calHi - calLo) * dt,
            gravityMag = gMag,
            accSigma = sigma,
            gyroBiasNorm = bias.norm(),
            dt = dt,
            gyroInstability = inst[0],
            gyroInstTauS = inst[1],
            attBudgetM = attBudget,
            vPeakMps = vPeak,
            anchorSuspect = anchorSuspect,
        )

        // ---- 12. 模型敏感 → 拒答（不再"给数 + 点名"）----
        // 原则是"宁可拒答，也不要一个像样的假数"。modelSensitive 的物理含义本来就是
        // "这个数由模型假设决定、不是由数据决定" —— 那它不该给数。
        // 原先照给再点名，实测后果：143734 报 1.44 cm 并附带一句"不可当真"，
        // 但 1.44 cm 是个**体面的**数字，用户没有真值可核对，只会照着读。
        // 全库 32 条给数里有 8 条落进这里，其中**只有 143734 有独立真值，而它是错的**；
        // 另一侧已知正确的 5 条（143719/143725/100444/100459/100529）一条都不在名单里。
        // 用 copy 而不是另写一份构造：拒答与给数的字段完全一致（含平台），
        // 被拒的那条自己仍然看得清，两边实现也好逐位对照。
        if (modelSensitive) {
            return out.copy(
                ok = false,
                failReason = String.format(Locale.US,
                    "这个数由残差模型决定、不是由数据决定（裸积分 %.2f cm、去趋势 %.2f cm，" +
                        "修正量是答案的 %.0f%%）—— 平台窄只说明端点稳，说明不了模型对。" +
                        "请重滑一次：快一点、起手和收手都停稳",
                    bare * 100, distHead * 100, modelGap / max(distHead, 1e-9) * 100),
            )
        }
        return out
    }
}

// ------------------------------------------------------------------ 静止检测

/** 逐样本局部标准差（居中窗，逐轴取最大）。居中 → 运动边界不滞后。 */
private fun rollingMaxStd(x: DoubleArray, y: DoubleArray, z: DoubleArray,
                          dt: Double, winS: Double): DoubleArray {
    val n = x.size
    val h = max(1, (winS / dt).toInt() / 2)
    val out = DoubleArray(n)
    for (i in 0 until n) {
        val a = max(0, i - h)
        val b = min(n, i + h + 1)
        val c = (b - a).toDouble()
        var sx = 0.0; var sy = 0.0; var sz = 0.0
        for (k in a until b) { sx += x[k]; sy += y[k]; sz += z[k] }
        val mx = sx / c; val my = sy / c; val mz = sz / c
        var vx = 0.0; var vy = 0.0; var vz = 0.0
        for (k in a until b) {
            val dx = x[k] - mx; vx += dx * dx
            val dy = y[k] - my; vy += dy * dy
            val dz = z[k] - mz; vz += dz * dz
        }
        out[i] = maxOf(sqrt(vx / c), sqrt(vy / c), sqrt(vz / c))
    }
    return out
}

/** 逐样本局部最大值（居中窗）。 */
private fun rollingMax(v: DoubleArray, dt: Double, winS: Double): DoubleArray {
    val n = v.size
    val h = max(1, (winS / dt).toInt() / 2)
    val out = DoubleArray(n)
    for (i in 0 until n) {
        val a = max(0, i - h)
        val b = min(n, i + h + 1)
        var m = 0.0
        for (k in a until b) if (v[k] > m) m = v[k]
        out[i] = m
    }
    return out
}

/** 与 Python `sorted(xs)[int(q*len)]` 一致。 */
private fun percentile(xs: DoubleArray, q: Double): Double {
    if (xs.isEmpty()) return 0.0
    val s = xs.clone()
    s.sort()
    val idx = min(s.size - 1, max(0, (q * s.size).toInt()))
    return s[idx]
}

/** 返回包含下标 [i] 的最大连续静止段 [start, end)。 */
private fun runContaining(flags: BooleanArray, i: Int, n: Int): IntArray {
    var a = i
    while (a > 0 && flags[a - 1]) a--
    var b = i + 1
    while (b < n && flags[b]) b++
    return intArrayOf(a, b)
}

/**
 * 静止段内估计陀螺零偏**不稳定性**（不是零偏本身）。
 *
 * 零偏的常值部分能被校准段平均掉；抹不掉的是它在测量时长内的游走，
 * 而那正是姿态漂移 → 重力泄漏的源头。做法：把静止段切成 4 段，
 * 相邻段均值之差的模长 /√2 就是该时长尺度上零偏变化的 σ。
 *
 * 返回 (σ rad/s, τ s)；τ 是估计所依据的时间尺度，外推到其他时长用 √(T/τ)。
 * 静止段太短（<16 样本）时返回 (0, 0) = 不可估。
 */
private fun gyroBiasInstability(
    gx: DoubleArray, gy: DoubleArray, gz: DoubleArray,
    run: IntArray, dt: Double,
): DoubleArray {
    val a = run[0]; val b = run[1]; val n = b - a
    if (n < 16) return doubleArrayOf(0.0, 0.0)
    val h = max(4, n / 4)
    val k = n / h
    if (k < 2) return doubleArrayOf(0.0, 0.0)
    val segs = Array(k) { j ->
        var sx = 0.0; var sy = 0.0; var sz = 0.0
        for (i in a + j * h until a + (j + 1) * h) {
            sx += gx[i]; sy += gy[i]; sz += gz[i]
        }
        V3(sx / h, sy / h, sz / h)
    }
    var acc = 0.0
    for (j in 0 until k - 1) acc += (segs[j + 1] - segs[j]).norm()
    return doubleArrayOf(acc / (k - 1) / Math.sqrt(2.0), h * dt)
}

/**
 * 姿态漂移造成的位移误差上界（米）。
 * 角速度误差 b_g 常值时：θ = b_g·t，泄漏 = g·θ，位移误差 = (1/6)·g·b_g·T³。
 * 零偏是游走的，在 T 尺度上的不确定度由 τ 尺度的 σ 按 √(T/τ) 扩散。
 */
private fun attitudeDriftBudget(gMag: Double, sigmaBg: Double, tau: Double, t: Double): Double {
    if (sigmaBg <= 0.0 || tau <= 0.0 || t <= 0.0) return 0.0
    val sigmaT = sigmaBg * Math.sqrt(t / tau)
    return gMag * sigmaT * t * t * t / 6.0
}

// ------------------------------------------------------------------ 姿态

/**
 * 互补滤波。**加速度修正在静止段才开**，运动中纯陀螺外推。
 *
 * v1 的参考是 `g_ref_body`（校准时刻的体重力方向，体坐标里固定不动）—— 等于让滤波器主动
 * 抵抗手机的真实旋转：持续角速率 ω 下稳态姿态误差 ≈ ω·dt/α（α=0.05、dt=9.5 ms 时
 * ≈0.19·ω），ω=5°/s 就有 1° 误差、重力泄漏 0.17 m/s²。
 * v2 换成当前**实测**加速度方向（标准 Mahony），手机转它就跟着转，没有这个稳态误差。
 *
 * [i0] 从第 i0 个样本才开始积分（之前的填初始姿态）。必须传校准段起点：初始姿态取自
 * 校准段的重力方向，若从样本 0 积起、而手机在校准段之前转过（实测 143734 转了 21°），
 * 这个初始姿态对样本 0 就是错的，而且一路带进窗口 —— 静止段残余达 3.3 m/s²。
 *
 * [gMag] 给了就用**模长**检验代替方向门限：方向门限会自锁 —— 误差一旦超过 GATE_RAD
 * 就永远"跳过修正"，再也拉不回来（143734 实测 19.5° 全程没修正过一次）。
 * 静止与否已由 flags 判定，模长检验表达的是同一件事的另一面：静止时 |acc| 就该 = |g|。
 */
private fun solveAttitude(
    ax: DoubleArray, ay: DoubleArray, az: DoubleArray,
    gx: DoubleArray, gy: DoubleArray, gz: DoubleArray,
    bias: V3, gRef: V3, flags: BooleanArray, dt: Double,
    i0: Int = 0, gMag: Double = 0.0,
): Array<Quat> {
    val n = ax.size
    val qs = Array(n + 1) { Quat.IDENTITY }
    var q = qFromTwoVectors(gRef, V3.UP)
    val start = i0.coerceIn(0, (n - 1).coerceAtLeast(0))
    for (k in 0..start) qs[k] = q       // qs[i] 对应 acc[i]：前 i0+1 个都是初始姿态
    for (i in start until n - 1) {
        val w = V3(gx[i] - bias.x, gy[i] - bias.y, gz[i] - bias.z)
        val qPred = q.mul(qDelta(w, dt)).normalized()
        if (flags[i]) {
            val e = qPred.rotateInv(V3.UP).cross(V3(ax[i], ay[i], az[i]).unit())
            val good = if (gMag > 0.0)
                abs(V3(ax[i], ay[i], az[i]).norm() - gMag) <= GATE_MAG_FRAC * gMag
            else e.norm() <= GATE_RAD
            q = if (!good) qPred
            else qPred.mul(qDelta(e * -ALPHA, dt)).normalized()
        } else {
            q = qPred
        }
        qs[i + 1] = q
    }
    return qs
}

// ------------------------------------------------------------------ 积分与 ZUPT

/** 末端位置随 Zupt 一起返回：杠杆臂修正是加在**位移矢量**上的，只有标量不够用。 */
private class Zupt(
    val distanceM: Double,
    val anchorDrift: Double,
    val endX: Double, val endY: Double, val endZ: Double,
)

/** 裸积分（不去趋势）：梯形二次积分，返回净位移模长。作为残差模型分歧的另一端。 */
/** 裸积分的末端**位移矢量**（不出模长）—— 杠杆臂修正必须加在矢量上，标量没法加。 */
private fun integrateEndVec(lx: DoubleArray, ly: DoubleArray, lz: DoubleArray,
                            w0: Int, w1: Int, dt: Double): DoubleArray {
    val m = w1 - w0
    var px = 0.0
    var py = 0.0
    var pz = 0.0
    var vx = 0.0
    var vy = 0.0
    var vz = 0.0
    for (k in 1..m) {
        val i = w0 + k
        val nx = vx + 0.5 * (lx[i - 1] + lx[i]) * dt
        val ny = vy + 0.5 * (ly[i - 1] + ly[i]) * dt
        val nz = vz + 0.5 * (lz[i - 1] + lz[i]) * dt
        px += 0.5 * (vx + nx) * dt
        py += 0.5 * (vy + ny) * dt
        pz += 0.5 * (vz + nz) * dt
        vx = nx; vy = ny; vz = nz
    }
    return doubleArrayOf(px, py, pz)
}

/**
 * 窗口内裸速度在**运动段** [m0, m1] 上的峰值（米/秒）。
 *
 * 只服务于 [MeasureResult.anchorSuspect]：如果标为「静止」的锚点段里那个速度，
 * 和运动段的速度峰值同量级，说明锚点段其实还在动。与 [integrateEndVec] 同一套梯形积分。
 */
private fun motionVelPeak(lx: DoubleArray, ly: DoubleArray, lz: DoubleArray,
                          w0: Int, w1: Int, m0: Int, m1: Int, dt: Double): Double {
    var vx = 0.0
    var vy = 0.0
    var vz = 0.0
    var pk = 0.0
    for (k in 1..(w1 - w0)) {
        val i = w0 + k
        vx += 0.5 * (lx[i - 1] + lx[i]) * dt
        vy += 0.5 * (ly[i - 1] + ly[i]) * dt
        vz += 0.5 * (lz[i - 1] + lz[i]) * dt
        if (i in m0..m1) {
            val s = sqrt(vx * vx + vy * vy + vz * vz)
            if (s > pk) pk = s
        }
    }
    return pk
}

private fun integrateEnd(lx: DoubleArray, ly: DoubleArray, lz: DoubleArray,
                         w0: Int, w1: Int, dt: Double): Double {
    val m = w1 - w0
    var px = 0.0
    var py = 0.0
    var pz = 0.0
    var vx = 0.0
    var vy = 0.0
    var vz = 0.0
    for (k in 1..m) {
        val i = w0 + k
        val nx = vx + 0.5 * (lx[i - 1] + lx[i]) * dt
        val ny = vy + 0.5 * (ly[i - 1] + ly[i]) * dt
        val nz = vz + 0.5 * (lz[i - 1] + lz[i]) * dt
        px += 0.5 * (vx + nx) * dt
        py += 0.5 * (vy + ny) * dt
        pz += 0.5 * (vz + nz) * dt
        vx = nx; vy = ny; vz = nz
    }
    return sqrt(px * px + py * py + pz * pz)
}

/**
 * 双锚点 ZUPT：把两侧静止段内**所有**样本的速度均值拉到 0。
 *
 * 设 v(t) = v_true(t) + A·t + B（残差加速度被当作常量）。A、B 由两个静止段的速度均值解出：
 * 取段内平均时刻 T 与平均速度 V，解 A·T+B = V。比只钉两个端点稳 √N，对单点毛刺也不敏感。
 *
 * 试过把残差模型加到二次（A·t + B·t²/2，允许残差加速度线性漂移）：仿真里偏差反而从
 * ±0.4% 涨到 −1~−2.6%，平台宽度涨 3~4 倍 —— 6 个参数用 6 个约束定，窗口一短矩阵就病态。
 * 所以只留线性。
 */
private fun zuptTwoAnchor(
    lx: DoubleArray, ly: DoubleArray, lz: DoubleArray,
    w0: Int, w1: Int, run0In: IntArray, run1In: IntArray, dt: Double,
): Zupt {
    val m = w1 - w0
    val vx = DoubleArray(m + 1); val vy = DoubleArray(m + 1); val vz = DoubleArray(m + 1)
    for (k in 1..m) {
        val i = w0 + k
        vx[k] = vx[k - 1] + 0.5 * (lx[i - 1] + lx[i]) * dt
        vy[k] = vy[k - 1] + 0.5 * (ly[i - 1] + ly[i]) * dt
        vz[k] = vz[k - 1] + 0.5 * (lz[i - 1] + lz[i]) * dt
    }
    // 锚点段裁到窗口内（平台扫描会传入比窗口更宽的段）
    var a0 = max(run0In[0], w0); var b0 = min(run0In[1], w1 + 1)
    var a1 = max(run1In[0], w0); var b1 = min(run1In[1], w1 + 1)
    if (b0 <= a0) { a0 = w0; b0 = w0 + 1 }
    if (b1 <= a1) { a1 = w1; b1 = w1 + 1 }

    fun meanT(a: Int, b: Int): Double {
        var s = 0.0
        for (i in a until b) s += (i - w0) * dt
        return s / (b - a)
    }
    val t0 = meanT(a0, b0)
    val t1 = meanT(a1, b1)
    val dt01 = t1 - t0
    val v0x = meanRange(vx, a0 - w0, b0 - w0)
    val v0y = meanRange(vy, a0 - w0, b0 - w0)
    val v0z = meanRange(vz, a0 - w0, b0 - w0)
    val v1x = meanRange(vx, a1 - w0, b1 - w0)
    val v1y = meanRange(vy, a1 - w0, b1 - w0)
    val v1z = meanRange(vz, a1 - w0, b1 - w0)
    val axErr = if (abs(dt01) > 1e-12) (v1x - v0x) / dt01 else 0.0
    val ayErr = if (abs(dt01) > 1e-12) (v1y - v0y) / dt01 else 0.0
    val azErr = if (abs(dt01) > 1e-12) (v1z - v0z) / dt01 else 0.0
    val bxErr = v0x - axErr * t0
    val byErr = v0y - ayErr * t0
    val bzErr = v0z - azErr * t0

    val px = DoubleArray(m + 1); val py = DoubleArray(m + 1); val pz = DoubleArray(m + 1)
    val cx = DoubleArray(m + 1); val cy = DoubleArray(m + 1); val cz = DoubleArray(m + 1)
    for (k in 0..m) {
        val t = k * dt
        cx[k] = vx[k] - (axErr * t + bxErr)
        cy[k] = vy[k] - (ayErr * t + byErr)
        cz[k] = vz[k] - (azErr * t + bzErr)
        if (k > 0) {
            px[k] = px[k - 1] + 0.5 * (cx[k - 1] + cx[k]) * dt
            py[k] = py[k - 1] + 0.5 * (cy[k - 1] + cy[k]) * dt
            pz[k] = pz[k - 1] + 0.5 * (cz[k - 1] + cz[k]) * dt
        }
    }
    val dist = sqrt(px[m] * px[m] + py[m] * py[m] + pz[m] * pz[m])
    // 锚点漂移：**修正前**速度在两侧静止段里的均值模长 —— 也就是"要去掉多少速度"。
    // （早先这里算的是修正**后**的均值，那是构造上恒等于 0 的死字段。）
    // 它和 |b|·T 同量级时才说明残差是常量；远大于 |b|·T 说明静止段里其实还在动。
    val drift = max(norm3(v0x, v0y, v0z), norm3(v1x, v1y, v1z))
    return Zupt(dist, drift, px[m], py[m], pz[m])
}

private fun meanRange(v: DoubleArray, from: Int, to: Int): Double {
    if (to <= from) return 0.0
    var s = 0.0
    for (i in from until to) s += v[i]
    return s / (to - from)
}

private fun norm3(x: Double, y: Double, z: Double) = sqrt(x * x + y * y + z * z)

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

// ------------------------------------------------------------------ 算法原语

private fun qDelta(w: V3, dt: Double) =
    Quat(1.0, w.x * dt / 2.0, w.y * dt / 2.0, w.z * dt / 2.0).normalized()

/** 求 q，使 R(q)·[aBody] == [bWorld]（两向量都当单位向量用）。 */
private fun qFromTwoVectors(aBody: V3, bWorld: V3): Quat {
    val a = aBody.unit()
    val b = bWorld.unit()
    val c = a.cross(b)
    val d = a.dot(b)
    if (d < -0.999999) {
        var axis = a.cross(V3(1.0, 0.0, 0.0))
        if (axis.norm() < 1e-6) axis = a.cross(V3(0.0, 1.0, 0.0))
        axis = axis.unit()
        return Quat(0.0, axis.x, axis.y, axis.z)
    }
    return Quat(1.0 + d, c.x, c.y, c.z).normalized()
}

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
private fun median(xs: DoubleArray): Double {
    if (xs.isEmpty()) return 0.0
    val s = xs.clone()
    s.sort()
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
