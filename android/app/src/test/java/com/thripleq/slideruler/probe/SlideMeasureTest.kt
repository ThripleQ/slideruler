package com.thripleq.slideruler.probe

import org.junit.Assert.assertEquals
import org.junit.Assert.assertFalse
import org.junit.Assert.assertTrue
import org.junit.Test

/**
 * 移植正确性：Kotlin 引擎必须逐位复现 Python `sim/measure2.py`。
 *
 * 为什么非要做这一层：v1 引擎曾经"看起来对" —— 两套实现都能跑出接近 15 cm 的数字，
 * 但判据差一点点换一条采集就崩（v1 真机上报过 168 cm、以及静止不动报 18.2 cm）。
 * 所以断言的不只是距离，还有**下标和分支选择**：窗口落在哪、锚点在哪、有没有拒答。
 * 那才是判据本身。
 *
 * 夹具是真机采集原件（acc.csv / gyr.csv / meta.txt），期望值由 Python 在同一份 CSV 上产出。
 * 覆盖七条互不相同的路径：
 *   100631  完全静止          → 判"没动"，距离 0
 *   100444  正常滑动（短停顿） → 15.4156 cm，锚点窄，裸积分与去趋势只差 3%
 *   100459  正常滑动（长停顿） → 15.7530 cm，锚点宽、平台较宽、修正量 46%（未超警戒线）
 *   100409  松手时还在动      → 6.8475 cm，修正量是答案的 2.6 倍 → **模型敏感**
 *   090006  尾部没静止段      → 必须**拒答**而不是给个像样的数（v1 在这里给过 15.91）
 *   100555  尾部静止只有 0.048 s → 必须**拒答**（旧版给 1.01 cm + 警告）
 *   100601  按下前一直在动    → 拒答，且理由必须指向"按下之前"而不是"松手之后"
 *   合成      按钮窗与运动零重叠  → 拒答（旧版会去量按下前挪手机那一段）
 *   合成      校准段被运动污染    → 拒答（这一条 Python 有、Kotlin 原本没有）
 *   122416  一次按住被切成 5 段  → 必须**合并**成 #256→#638（旧版只取 0.29 s，报 0.59 cm）
 *   122303  窗内两段运动        → 合并成 #176→#503，平台仍窄（证明合并不是一律搞炸）
 *
 * 下标约定：`measure2.py` 内部把连续段存成半开区间 `[a,b)`，但它另外导出一个
 * `anchorInclusive` = `((a0,b0-1),(a1,b1-1))`。Kotlin 侧统一报**闭区间**（与
 * `motionIndex`/`windowStart`/`windowEnd` 一致），所以这里比对的是 `anchorInclusive`。
 *
 * ⚠ 这个 1e-6 m 容差检验的是「Kotlin 复现了 Python 的同一个答案」，
 * **不是**「测量误差小于百万分之一米」。后者需要硬真值，见 docs/v2-comparison.md 的残留一节。
 */
class SlideMeasureTest {

    private class Fixture(val acc: List<Pt3>, val gyr: List<Pt3>, val downNs: Long, val upNs: Long)

    private fun load(name: String): Fixture {
        val dir = "fixtures/$name"
        val down = Regex("DOWN\\s+event=\\d+ ms\\s+recv_uptime=\\d+ ms\\s+recv_elapsed=(\\d+) ns")
            .find(read(dir, "meta.txt"))!!.groupValues[1].toLong()
        val up = Regex("UP\\s+event=\\d+ ms\\s+recv_uptime=\\d+ ms\\s+recv_elapsed=(\\d+) ns")
            .find(read(dir, "meta.txt"))!!.groupValues[1].toLong()
        return Fixture(readCsv(dir, "acc.csv"), readCsv(dir, "gyr.csv"), down, up)
    }

    private fun read(dir: String, file: String): String =
        javaClass.classLoader!!.getResourceAsStream("$dir/$file")!!.bufferedReader().readText()

    private fun readCsv(dir: String, file: String): List<Pt3> =
        read(dir, file).lineSequence().drop(1)
            .filter { it.isNotBlank() }
            .map {
                val p = it.trim().split(",")
                Pt3(p[0].toLong(), p[1].toDouble(), p[2].toDouble(), p[3].toDouble())
            }.toList()

    // 1 m 的 1e-6 = 0.0001 cm —— 比界面显示的 0.01 cm 还紧 100 倍
    private val tol = 1e-6

    /**
     * 测量头几何（杠杆臂）：Δp_H = Δp_I + [R(t1)−R(t0)]·r
     *
     * 期望值由 `sim/measure2.py` 在同一份夹具上以 head_r=(0,−0.09,0) 产出。
     * 这条断言锁的是**两边对同一个几何修正的理解一致**，不是"加了杠杆臂就更准" ——
     * r 必须标定，见 docs/v4-head-model.md §3。
     */
    @Test
    fun `杠杆臂 - 测量头位移 = IMU 位移 + 两端姿态差乘 r`() {
        val f = load("20261011_100444")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs, 0.0, -0.09, 0.0)
        assertTrue(r.ok)
        assertEquals("测量头位移", 0.14426168222690888, r.distanceM, tol)
        assertEquals("IMU 位移（应与不带 headR 时一致）", 0.15415506077049201, r.imuDistM, tol)
        assertEquals("杠杆臂项", 0.0099533247435707912, r.leverM, tol)
        assertEquals("裸积分也要加杠杆臂项", 0.14278939325128673, r.bareM, 1e-6)
        assertEquals("平台下界（端点在动，杠杆臂一起扫）", 0.14385469818057917, r.plateauLoM, tol)
        assertEquals("平台上界", 0.14426211040095777, r.plateauHiM, tol)

        // 关掉杠杆臂时必须逐位回到 v3 的答案 —— 这是"默认不改变行为"的守卫
        val v3 = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertEquals("headR=0 时 distanceM 必须等于 imuDistM", v3.distanceM, v3.imuDistM, 0.0)
        assertEquals("headR=0 时杠杆臂项必须为 0", 0.0, v3.leverM, 0.0)
        assertEquals("headR=0 时必须复现 v3 的距离", 0.15415506077049201, v3.distanceM, tol)
    }

    /**
     * 起手前静止段的**角抖动** —— 它按场次分层，但**不是**陀螺零偏游走。
     *
     * 这条是 2026-10-16 那次证伪的守卫。库里 `20261011_085650` 是一段 76 s 的纯静止
     * 录制（meta 按住时长 0 ms，没做滑动），在那上面：
     *   · 5 s 窗口的零偏均值全部落在 ±3e-5 rad/s 内，随机散布、无单调漂移；
     *   · Allan deviation 从 10 ms 到 9.7 s **一路单调下降**（0.00052 → 7.8e-6），
     *     没有平台（flicker 的斜率 0）也没有回升（速率随机游走的 +0.5）。
     * ⟹ 这台 K40 的零偏**不游走**。而本字段在真机上给 0.0001~0.17、且按场次分层：
     *   · 桌面受控场次 100444 / 100529 → 0.00036 / 0.00012（≈ 器件噪声底）；
     *   · "姿态随意"手持场次 143719 / 143725 → 0.0123 / 0.0149（差 30~120 倍）。
     * ⟹ 它量的是"起手前手抖不抖"（端点参考稳不稳），不是器件漂移。
     *
     * 因此由它派生的拒答与报警**已全部删除**：分层只划开了两个场次，划不开同一
     * 场次里的对与错（143532 抖动 0.0067 给出 15.92 cm、143604 抖动 0.0021 给出
     * 14.42 cm，两条都准）。没有标注数据就定阈值 = 拿真值调参。
     * ⚠ 别再往这条上加 ∝T³ 的"姿态漂移上界"：它算出的预算与实际误差**反相关**
     * （预算 5.4 cm 的 143734 错得最狠 −90%，预算 35.7 cm 的 143725 反而 +5%）。
     * 期望值由 sim/measure2.py 在同一份夹具上产出。
     */
    @Test
    fun `起手前角抖动 - 按场次分层但不得派生出拒答或T3报警`() {
        val f = load("20261011_100444")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertTrue(r.ok)
        assertEquals("桌面场次的角抖动（≈器件噪声底）", 0.00035816660652023468, r.startTremorRadS, 1e-9)
        assertEquals("τ", 0.076109999999999997, r.startTremorTauS, 1e-9)

        // 手持场次：同一个量高出 34 倍，但读数照给 —— 抖动大不等于一定错
        val g = load("20261011_143719")
        val rg = SlideMeasure.measure(g.acc, g.gyr, g.downNs, g.upNs)
        assertTrue("抖动大不拒答 —— 判据没有干净分界线", rg.ok)
        assertEquals("手持场次的角抖动", 0.012304473000221985, rg.startTremorRadS, 1e-9)
        assertEquals("τ", 0.24736696400000002, rg.startTremorTauS, 1e-9)
        assertTrue("两个场次相差 >30 倍，是本字段唯一可断言的性质",
            rg.startTremorRadS > 30.0 * r.startTremorRadS)

        val h = load("20261011_100529")
        val rh = SlideMeasure.measure(h.acc, h.gyr, h.downNs, h.upNs)
        assertEquals("同场次里最低的一条", 0.00011947117157198642, rh.startTremorRadS, 1e-9)
    }

    /**
     * 按钮窗内禁止加速度修正 —— 别让「匀速段」污染姿态。
     *
     * 122416 是一次病态按住：静止检测把一次滑动切成了 5 段，而被判成静止的那些
     * **匀速段**样本，其 |‖acc‖−g| ≈ 0，**能通过模长门控**。
     * 姿态解算若只按 `flags` 放行（Python 用 `gate_flags` 额外排除按钮窗），
     * 它们就会被拿去"修正"姿态，把运动加速度写进重力方向。
     *
     * ⚠ 这是 v4 的关键改动，Kotlin 移植时漏了：修之前两侧在这条夹具上差 9.4 cm
     * （Kotlin 1.4688 / Python 1.5632），而**其他夹具全都看不出来** ——
     * 正常数据的运动段 |‖acc‖−g| 很大，门控本来就会拒，放行与不放行结果相同。
     * 所以这条夹具是这套判据唯一的守卫，期望值由 sim/measure2.py 产出。
     */
    @Test
    fun `按钮窗内禁止加速度修正 - 匀速段不许污染姿态`() {
        val f = load("20261011_122416")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertEquals("裸积分（姿态一旦被污染，这里最先变）",
            5.570952399021404, r.bareM, 1e-9)
        assertEquals("锚点漂移", 3.6199266330942566, r.anchorDrift, 1e-9)
        assertEquals("距离（只作诊断，不展示给用户）", 1.563247657506976, r.distanceM, 1e-6)
        assertEquals("区间选择与姿态无关，不受影响", 256, r.motionIndex)
        assertFalse("这个数由模型决定 → 拒答", r.ok)
    }

    @Test
    fun `完全静止 - 判没动而不是凭空造距离`() {
        val f = load("20261011_100631")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertTrue("应当成功（只是没动）", r.ok)
        assertTrue("应当判为「没有运动」", r.noMotion)
        assertFalse("不该是慢动作分支", r.slowMotion)
        assertEquals("距离必须是 0", 0.0, r.distanceM, 0.0)
        // v1 在这条采集上给 18.20 cm —— 这就是用户报的"按着不动也能测出很长一段距离"
    }

    @Test
    fun `正常滑动 - 距离与平台`() {
        val f = load("20261011_100444")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertTrue(r.ok)
        assertFalse(r.noMotion)
        assertEquals("距离", 0.15415506077049201, r.distanceM, tol)
        assertEquals("平台下界", 0.15374307803284132, r.plateauLoM, tol)
        assertEquals("平台上界", 0.15415500053000761, r.plateauHiM, tol)
        // 判据本身：运动段、窗口、两侧锚点
        assertEquals("运动段起点 #431", 431, r.motionIndex)
        assertEquals("窗口起点 #416", 416, r.windowStart)
        assertEquals("窗口终点 #676", 676, r.windowEnd)
        assertEquals("左锚点起点 #396", 396, r.anchorStart)
        assertEquals("右锚点终点 #689（闭区间）", 689, r.anchorEnd)
        assertEquals("运动时长", 2.1881625, r.motionS, 1e-6)
        assertEquals("窗口时长", 2.473575, r.windowS, 1e-6)
        assertEquals("校准段", 0.33298125, r.calS, 1e-6)
        assertEquals("松手后留下的静止", 0.27589874999999997, r.tailStaticS, 1e-6)
        assertFalse("尾部静止足够，不该报终点锚点问题", r.endHitBound)
        // 裸积分与去趋势只差 1% —— 残差模型分歧小，这个数由数据决定
        assertEquals("裸积分", 0.15259275303446607, r.bareM, 1e-6)
        assertEquals("模型分歧", 0.0015623077360259474, r.modelGapM, 1e-6)
        assertFalse("修正量远小于答案，不该报模型敏感", r.modelSensitive)
        assertEquals("锚点漂移", 0.01714653384141511, r.anchorDrift, 1e-9)
        assertEquals("实测重力模长", 9.913901720863374, r.gravityMag, 1e-9)
        assertEquals("静止段噪声地板", 0.013673987754797903, r.accSigma, 1e-9)
        assertTrue("平台应当很窄", r.confidence < 0.02)
    }

    @Test
    fun `长停顿滑动 - 停顿必须被裁掉`() {
        val f = load("20261011_100459")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertTrue(r.ok)
        assertEquals("距离", 0.15753443241336959, r.distanceM, tol)
        assertEquals("平台下界", 0.15688152894988547, r.plateauLoM, tol)
        assertEquals("平台上界", 0.16440202078475305, r.plateauHiM, tol)
        assertEquals(439, r.motionIndex)
        assertEquals(424, r.windowStart)
        assertEquals(723, r.windowEnd)
        assertEquals("左锚点其实是运动前那一大段静止 #275", 275, r.anchorStart)
        assertEquals("右锚点终点 #805（闭区间）", 805, r.anchorEnd)
        // 用户按住 1.87 s，但真实运动 2.56 s（松手后还在动）→ 窗口 2.84 s。
        // 关键是窗口**没有被拉成**运动前那 1.56 s 静止 —— 那段静止对距离零贡献、
        // 对误差是完整的 ∝T² 加成。
        assertEquals("按住时长", 1.87420875, r.holdS, 1e-6)
        assertEquals("真实运动时长", 2.5591987499999997, r.motionS, 1e-6)
        assertEquals("窗口时长必须接近运动时长", 2.84461125, r.windowS, 1e-6)
        assertEquals("校准段取的是停顿那一段", 1.560255, r.calS, 1e-6)
        assertEquals("松手后留下的静止", 0.9323475, r.tailStaticS, 1e-6)
        assertTrue("平台偏宽 → 可信度应当一般", r.confidence > 0.02)
    }

    @Test
    fun `平台宽度是自带的误差棒`() {
        val f = load("20261011_100409")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertTrue(r.ok)
        assertEquals("距离", 0.068496449074256685, r.distanceM, tol)
        assertEquals("平台下界", 0.063128874735805085, r.plateauLoM, tol)
        assertEquals("平台上界", 0.076194522295233283, r.plateauHiM, tol)
        assertEquals(400, r.motionIndex)
        assertEquals(385, r.windowStart)
        assertEquals(694, r.windowEnd)
        assertEquals(353, r.anchorStart)
        assertEquals(707, r.anchorEnd)
        assertEquals("松手后留下的静止", 0.27589724199999999, r.tailStaticS, 1e-6)
        // 平台 (7.61-6.31)/6.85 = 19% —— 换个端点答案就变，这数本来就不该太当真
        assertTrue("平台占比应当很大", r.confidence > 0.15)
    }

    /**
     * 修正量超过答案本身 —— 这个数由**模型假设**决定，平台窄救不了它。
     *
     * 122416 的裸积分 537 cm 被去趋势修成 147 cm，修正量是答案的 2.7 倍：
     * 答案是两个大数之差，换一个残差模型就换一个答案。
     *
     * ⚠ 这条以前用的是 100409。修掉姿态初始化错位（校准段之前转过的角度没被计入）
     * 之后，100409 的修正量从答案的 260% 降到 26%，不再触发 —— 也就是说那条数据
     * 原先的"模型敏感"其实是姿态 bug 的症状，不是模型分歧。判据本身没变，
     * 是它该报的那批数据变少了：真机 32 条里从 20 条降到 8 条。
     *
     * 立这条判据的依据是真机分界很干净：三条有硬真值的是 1% / 4% / 3% 都没超，
     * 而没有真值的那批（150% / 260% / 590% / 4587%）全超。见 sim/compare_v1_v2.py。
     */
    @Test
    fun `修正量超过答案本身 - 必须拒答而不是给一个体面的数`() {
        val f = load("20261011_122416")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        // 原先这里是 assertTrue(r.ok)，然后给出 147 cm 再挂一句"不可当真"。
        // 现在改成拒答：点名的同时给数，实测后果是用户照着读 —— 143734 报 1.44 cm
        // 就是个**体面的**错数（1.44 看起来完全正常，而真值是 15）。
        assertFalse("修正量是答案的 2.7 倍 → 必须拒答", r.ok)
        assertTrue("拒答也要由模型敏感旗标点名", r.modelSensitive)
        assertTrue(
            "理由要说清是模型决定的：${r.failReason}",
            r.failReason!!.contains("残差模型决定"),
        )
        // 拒答了也要把证据留下，否则这条数据就再也看不清楚了
        assertEquals("裸积分", 5.570952399021404, r.bareM, 1e-6)
        assertEquals("模型分歧", 4.007704741514427, r.modelGapM, 1e-6)
        assertEquals("锚点漂移", 3.6199266330942566, r.anchorDrift, 1e-9)
        assertEquals("运动段速度峰值", 3.4154004650754186, r.vPeakMps, 1e-9)
        assertTrue("3.62 vs 3.42 m/s —— 锚点也在动", r.anchorSuspect)
        assertEquals("平台照算（拒答的条同样看得清）",
            1.2122176967837292, r.plateauLoM, 1e-6)

        // 对照：修正量只有 21% 的那条必须继续给数 —— 判据不能把所有去趋势都当成可疑
        val g = load("20261011_100459")
        val rg = SlideMeasure.measure(g.acc, g.gyr, g.downNs, g.upNs)
        assertTrue("100459 必须给数", rg.ok)
        assertFalse("21% < 100%，不该报模型敏感", rg.modelSensitive)
        assertFalse("锚点没在动（0.025 vs 0.155 m/s）", rg.anchorSuspect)
        assertEquals("裸积分", 0.12411948074406344, rg.bareM, 1e-6)
        assertEquals("模型分歧", 0.03341495166930615, rg.modelGapM, 1e-6)
    }

    /**
     * 终点锚点退化 —— 现在**必须拒答**，不能给数。
     *
     * 100555 里运动一直延续到 #971，而数据只有 976 个样本 —— 留给终点锚点的静止只有
     * 0.048 s，锚点被裁到 1~4 个样本。以前它会给 1.01 cm（平台 0.90~1.51，宽 60%）
     * 加一句警告；但带警告的距离，用户照样会当成答案。
     *
     * 为什么是"整段静止"而不是"裁进窗口的余量"：后者恒等于 MARGIN_S（0.143 s / 15 样），
     * 拿它去卡 STATIC_MIN_S=0.20 会把全部测量干掉。而想让锚点变长只能把窗口往外扩 ——
     * sim/exp_anchor.py 试过，三条有真值的里两条变差（+5.02%→+9.37%、+3.72%→+6.04%），
     * 偏差随窗口长度单调增长。所以锚点就这么长，这里只挡"假静止"。
     */
    @Test
    fun `尾部静止不够 - 必须拒答而不是给个带警告的数`() {
        val f = load("20261011_100555")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertFalse("静止段只有 0.04 s，不配当锚点，必须拒答", r.ok)
        val why = r.failReason!!
        assertTrue("要说清是松手之后那一侧：$why", why.contains("松手之后"))
        assertTrue("要说出实际静止时长：$why", why.contains("0.04"))
        assertTrue("要说出下限：$why", why.contains("0.20"))
    }

    /**
     * 左侧缺锚点时的**理由必须指对方向**。
     *
     * 100601 的运动从 #34 一路延续到 #339（正好是搜索区起点 → 按下前 2 s 手机一直在动），
     * 静止段只剩松手后那一段。用户再怎么「松手后保持不动」也修不好这个 ——
     * 他得先放稳再按。所以提示不能写「松手后没停稳」，那会把人引到错的地方去。
     */
    @Test
    fun `按下前手机一直在动 - 拒答理由必须指向正确的一侧`() {
        val f = load("20261011_100601")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertFalse("左侧没有静止段，必须拒答", r.ok)
        val why = r.failReason!!
        assertTrue("要指出是按下之前的问题：$why", why.contains("按下之前"))
        assertFalse("不能把锅甩给松手：$why", why.contains("没停稳"))
    }

    @Test
    fun `尾部没有静止段 - 必须拒答`() {
        val f = load("20261011_090006")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertFalse("取不到两端静止就必须拒答", r.ok)
        assertTrue(
            "理由要说清是「松手后没停稳」：${r.failReason}",
            r.failReason!!.contains("没停稳"),
        )
    }

    /**
     * 「取与按钮窗重叠最多的运动段」在**所有段都不重叠**时会退化成「取时间最早的一段」。
     *
     * 100444 的搜索区里有三段运动：#90-212、#294-395、#431-642。把按钮窗放到
     * #220~#280（两段运动之间的空隙），三段与它的重叠**全是 0**，
     * 旧代码于是选中了最早那一段（用户按下之前挪手机），量出 0.06 cm。
     *
     * 注意 `UP <= DOWN` 那个守卫挡不住这里：它只保证按钮窗长度为正，
     * 保证不了选中的运动就是这次滑动。所以必须单独查重叠 > 0。
     */
    @Test
    fun `按钮窗与所有运动段都不重叠 - 必须拒答`() {
        val f = load("20261011_100444")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.acc[220].tNs, f.acc[280].tNs)
        assertFalse("按钮期间没有运动，必须拒答", r.ok)
        assertTrue(
            "理由要说清是按下期间没测到运动：${r.failReason}",
            r.failReason!!.contains("按下期间没有检测到运动"),
        )
    }

    /**
     * 校准段混进了运动样本 —— Python 会拒，Kotlin 原本不会（移植时漏掉的一条）。
     *
     * 混着运动的均值既不是静止零偏也不是运动加速度，拿它当重力方向会让整条链路
     * 漏掉一大块重力。实测 100339：校准段被污染后 |g| 读到 9.43（本机应 9.91）、
     * σ=2.01，v1 据此报了 168.62 cm。
     *
     * 这一条之所以长期没被发现，是因为夹具里没有覆盖"校准污染"这个分支 ——
     * 单测的 1e-6 容差只能验证两边算出**同一个答案**，验证不了两边**走同一条分支**。
     */
    @Test
    fun `校准段被运动污染 - 必须拒答`() {
        val f = load("20261011_100444")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.acc[650].tNs, f.acc[700].tNs)
        assertFalse("校准段静止占比只有 80%，必须拒答", r.ok)
        val why = r.failReason!!
        assertTrue(
            "理由要说清是校准段不够静止：$why",
            why.contains("校准段不够静止") && why.contains("80"),
        )
    }

    /**
     * 一次按住里的多段运动必须合并，不许把一次滑动切碎。
     *
     * 122416 是一次真实按住（2.81 s）。静止检测在窗内切出 5 段运动
     * （0.20 / 0.29 / 0.14 / 0.26 / 1.19 s），旧逻辑"取与按钮窗重叠最多的一组"
     * 挑中了 0.29 s 那一段 → 报 **0.59 cm**。为什么挑错：真正那次 1.19 s 的滑动
     * 延续到松手之后，它在按钮窗内的重叠只有 0.24 s，输给了完全落在窗内的 0.29 s 抖动。
     *
     * 这里锁的是**区间选择**，不是距离精度：合并后必须覆盖 #256→#638。
     * 至于这个数本身 —— 平台 121~323 cm、修正量是答案的 2.6 倍，
     * `modelSensitive` 必须为真，所以这条现在**拒答**（只留诊断用的距离）。
     * 报一个明摆着离谱的数，好过报一个 0.59 cm 这种"看着体面但错"的数。
     */
    @Test
    fun `一次按住里的多段运动必须合并 - 不许切碎滑动`() {
        val f = load("20261011_122416")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        // 区间选择是这条测试真正锁的东西 —— 拒答不影响它，字段照旧带着
        assertEquals("合并后的运动区间起点", 256, r.motionIndex)
        assertEquals("窗口起点", 241, r.windowStart)
        assertEquals("窗口终点", 653, r.windowEnd)
        assertEquals("这次按住被静止检测切成了 5 段，必须全部合并", 5, r.mergedRuns)
        assertEquals("距离（只作诊断，不展示给用户）", 1.563247657506976, r.distanceM, tol)
        assertFalse("这个数由模型决定 → 拒答", r.ok)
        assertTrue("拒答由模型敏感旗标点名", r.modelSensitive)
    }

    /**
     * 终点锚点在动 —— 「松手」不等于「停住」。
     *
     * 143734（15 cm 硬真值）是本项目最贵的一条：**裸积分 16.19 cm 本来是对的**，
     * 双锚点扣掉 15.42 cm，报 1.44 cm。扣掉的不是漂移 —— 终点锚点段（标为「静止」的
     * 0.143 s）里积分速度稳定在 −0.31 m/s，而运动段速度峰值才 0.314 m/s。也就是
     * 手离开后手机还在匀速滑，而**匀速移动与静止在加速度计上不可分辨**（测的是比力，
     * 匀速时比力就是重力），静止检测必然把它判成静止。同一条的起点锚点是真静止
     * （|V0| = 0.007 m/s），问题只在终点一侧。
     *
     * 这条现在同时被两个判据抓到：模型敏感（→ 拒答）与锚点可疑（→ 报警）。
     * 锚点可疑**不拒答**（143725 同样可疑却给出 +5% 的正确定数，判据没有干净分界线，
     * 拿它拒答等于用真值调参），所以它必须能单独被观察到 —— 就是这条测试。
     */
    @Test
    fun `终点锚点在动 - 松手后手机还在匀速滑必须报出来`() {
        val f = load("20261011_143734")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertTrue("标为「静止」的终点锚点里速度与运动峰值同量级 → 必须点出来", r.anchorSuspect)
        assertEquals("终点锚点段里的速度", 0.3160358900303602, r.anchorDrift, 1e-9)
        assertEquals("运动段速度峰值", 0.3136389804577817, r.vPeakMps, 1e-9)
        assertTrue("两者同量级（比值 1.01）正是可疑的定义", r.anchorDrift > 0.5 * r.vPeakMps)
        assertFalse("这条同时被模型敏感抓到，所以拒答", r.ok)

        // 对照：同一次实验里正常的那三条，锚点都不该可疑。
        // 100444/100459/100529 的比值是 0.09 / 0.16 / 0.06，与 143 系列（0.90~1.01）分得很开。
        for (nm in listOf("20261011_100444", "20261011_100459", "20261011_100529")) {
            val g = load(nm)
            val rg = SlideMeasure.measure(g.acc, g.gyr, g.downNs, g.upNs)
            assertTrue("$nm 必须给数", rg.ok)
            assertFalse(
                "$nm 的终点锚点不该可疑（速度 ${rg.anchorDrift} vs 峰值 ${rg.vPeakMps}）",
                rg.anchorSuspect,
            )
        }
    }

    /**
     * 合并的温和情形：窗内有两段运动（0.20 s 的小抖动 + 2.41 s 的滑动），
     * 合并后 11.54 cm，平台 11.36~11.79 cm（很窄）→ 不报模型敏感。
     * 这条用来证明"合并"不是一律把结果搞炸：数据本身干净时它照样给出稳定的数。
     */
    @Test
    fun `合并两段运动 - 干净数据下平台仍然要窄`() {
        val f = load("20261011_122303")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertTrue(r.ok)
        assertEquals("合并后起点要包含那段小抖动", 176, r.motionIndex)
        assertEquals("窗内两段被合并", 2, r.mergedRuns)
        assertEquals("距离", 0.11536100593437933, r.distanceM, tol)
        // 注意：Python 的 plateau 导出的单位是 cm，Kotlin 的 plateauLoM/HiM 是 m
        assertEquals("平台下界 11.36cm", 0.11355627291419963, r.plateauLoM, tol)
        assertEquals("平台上界 11.79cm", 0.11790594818776837, r.plateauHiM, tol)
        assertFalse("干净数据不该报模型敏感", r.modelSensitive)
    }

    @Test
    fun `退化输入必须被拒绝而不是给个数`() {
        val f = load("20261011_100444")
        // DOWN==UP：所有运动组与按钮窗的重叠都是 0，"取重叠最大的一组"退化成"取最早的一组"，
        // 会去测按住之前那一段。Python 在这里报 0.06 cm（真的滑动是 15.4 cm）—— 答案看着无害，
        // 但它答的是另一个问题。必须拒。
        val same = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.downNs)
        assertFalse("DOWN == UP 必须拒答", same.ok)
        assertTrue("理由要说清是时刻异常：${same.failReason}", same.failReason!!.contains("时刻异常"))
        // UP 早于 DOWN 同理
        assertFalse(
            "UP < DOWN 必须拒答",
            SlideMeasure.measure(f.acc, f.gyr, f.upNs, f.downNs).ok,
        )
        assertFalse(
            "样本太少",
            SlideMeasure.measure(f.acc.take(10), f.gyr.take(10), f.downNs, f.upNs).ok,
        )
    }
}
