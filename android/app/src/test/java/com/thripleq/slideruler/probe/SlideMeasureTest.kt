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
 *   100555  尾部静止只有 0.048 s → 仍给数，但必须点名终点锚点退化、平台 60%
 *   100601  按下前一直在动    → 拒答，且理由必须指向"按下之前"而不是"松手之后"
 *
 * 下标约定：`measure2.py` 内部把连续段存成半开区间 `[a,b)`，但它另外导出一个
 * `anchorInclusive` = `((a0,b0-1),(a1,b1-1))`。Kotlin 侧统一报**闭区间**（与
 * `motionIndex`/`windowStart`/`windowEnd` 一致），所以这里比对的是 `anchorInclusive`。
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
        assertEquals("距离", 0.15415568278434755, r.distanceM, tol)
        assertEquals("平台下界", 0.15374404190282473, r.plateauLoM, tol)
        assertEquals("平台上界", 0.15415567595495647, r.plateauHiM, tol)
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
        // 裸积分与去趋势只差 3% —— 残差模型分歧小，这个数由数据决定
        assertEquals("裸积分", 0.1588530524781836, r.bareM, 1e-6)
        assertEquals("模型分歧", 0.004697369693836051, r.modelGapM, 1e-6)
        assertFalse("修正量远小于答案，不该报模型敏感", r.modelSensitive)
        assertEquals("锚点漂移", 0.023303692808520592, r.anchorDrift, 1e-9)
        assertEquals("实测重力模长", 9.913901720863374, r.gravityMag, 1e-9)
        assertEquals("静止段噪声地板", 0.013673987754797903, r.accSigma, 1e-9)
        assertTrue("平台应当很窄", r.confidence < 0.02)
    }

    @Test
    fun `长停顿滑动 - 停顿必须被裁掉`() {
        val f = load("20261011_100459")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertTrue(r.ok)
        assertEquals("距离", 0.15753003184370357, r.distanceM, tol)
        assertEquals("平台下界", 0.1568815289498855, r.plateauLoM, tol)
        assertEquals("平台上界", 0.1640571877892311, r.plateauHiM, tol)
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
        assertEquals("距离", 0.06847493525602662, r.distanceM, tol)
        assertEquals("平台下界", 0.06312887473580532, r.plateauLoM, tol)
        assertEquals("平台上界", 0.07605511186398839, r.plateauHiM, tol)
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
     * 100409 的裸积分 24.67 cm 被去趋势修成 6.85 cm，修正量是答案的 2.6 倍：
     * 答案是两个大数之差。此时平台扫描给出的 6.31~7.61（19%）看着还行，
     * 但平台只反映**端点抖动**，反映不了「残差到底是不是线性的」这个模型分歧。
     *
     * 立这条判据的依据是真机分界很干净：三条有真值（15 cm）的是 3% / 16% / 46% 都没超，
     * 四条没有真值的是 150% / 260% / 590% / 4587% 全超。见 sim/compare_v1_v2.py。
     * 它的误报率在仿真里极低（干净静止段 + 完美常量零偏，B/C/D/E 组 0/12）。
     */
    @Test
    fun `修正量超过答案本身 - 必须判为模型敏感而不是高可信`() {
        val f = load("20261011_100409")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertTrue(r.ok)
        assertEquals("裸积分", 0.2467165805266843, r.bareM, 1e-6)
        assertEquals("模型分歧", 0.17824164527065767, r.modelGapM, 1e-6)
        assertTrue("修正量是答案的 2.6 倍 → 必须报模型敏感", r.modelSensitive)
        assertEquals("锚点漂移", 0.2113703912072919, r.anchorDrift, 1e-9)

        // 对照：修正量 46%（< 100%）的那条不该触发 —— 判据不能把所有去趋势都当成可疑
        val g = load("20261011_100459")
        val rg = SlideMeasure.measure(g.acc, g.gyr, g.downNs, g.upNs)
        assertEquals("裸积分", 0.08527852130552778, rg.bareM, 1e-6)
        assertEquals("模型分歧", 0.072251510538175787, rg.modelGapM, 1e-6)
        assertFalse("46% < 100%，不该报模型敏感", rg.modelSensitive)
    }

    /**
     * 终点锚点退化：数据在手机刚停下时就到头了。
     *
     * 100555 里运动一直延续到 #971，而数据只有 976 个样本 —— 留给终点锚点的静止只有
     * 0.0476 s（< 余量 0.15 s），锚点被裁到 1 个样本。它给出的 1.01 cm 完全不值得信。
     * 以前这种情况只由平台宽度（1.51-0.90，60%）间接暴露；现在必须直接点名，
     * 而且要能区分「手机没停稳」和「录太短了」这两件事。
     */
    @Test
    fun `尾部静止不够 - 必须明说这个数不可靠`() {
        val f = load("20261011_100555")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)
        assertTrue("还是会给个数（但要把话说清）", r.ok)
        assertEquals("距离", 0.01009610388640131, r.distanceM, tol)
        assertEquals("留给终点锚点的静止", 0.04756849, r.tailStaticS, 1e-6)
        assertTrue("必须报出终点锚点有问题", r.endHitBound)
        assertTrue("尾部静止小于余量 0.15 s", r.tailStaticS < SlideMeasure.TAIL_MIN_SECONDS)
        assertTrue("平台应当极宽", r.confidence > 0.5)
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
