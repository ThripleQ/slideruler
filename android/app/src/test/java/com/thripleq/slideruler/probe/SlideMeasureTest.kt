package com.thripleq.slideruler.probe

import org.junit.Assert.assertEquals
import org.junit.Assert.assertFalse
import org.junit.Assert.assertTrue
import org.junit.Test

/**
 * Kotlin 测距引擎 vs Python 参考实现的逐位对齐。
 *
 * 移植最大的风险是「看起来对」：两套实现都能跑出接近 15 cm 的数字，但分支判据
 * （用不用去趋势、要不要顺延终点）差一点点就会在别的采集上崩掉。所以这里拿两条
 * **真实采集**做夹具，断言的不只是距离，还有下标和分支选择 —— 那些才是判据本身。
 *
 * 期望值由 `sim/slide_pipeline.py` 在同一份 CSV 上产出（见 docs/porting.md）。
 * 输入是同一份文件、同样的解析方式，所以两边应该几乎逐位相同；容差 1e-5 m
 * 只用来看住浮点求和顺序的差异，比显示的 0.01 cm 精度还紧 10 倍。
 */
class SlideMeasureTest {

    private class Fixture(
        val acc: List<Pt3>,
        val gyr: List<Pt3>,
        val downNs: Long,
        val upNs: Long,
    )

    private fun text(path: String): String =
        javaClass.getResourceAsStream(path)?.readBytes()?.toString(Charsets.UTF_8)
            ?: error("夹具缺失：$path")

    private fun readCsv(path: String): List<Pt3> {
        val out = ArrayList<Pt3>(2048)
        var header = true
        for (line in text(path).lineSequence()) {
            if (header) { header = false; continue }
            if (line.isBlank()) continue
            val f = line.split(',')
            if (f.size < 4) continue
            out.add(Pt3(f[0].trim().toLong(), f[1].trim().toDouble(), f[2].trim().toDouble(), f[3].trim().toDouble()))
        }
        return out
    }

    private fun load(name: String): Fixture {
        val meta = text("/fixtures/$name/meta.txt")
        val down = Regex("DOWN.*?recv_elapsed=(\\d+) ns").find(meta)!!.groupValues[1].toLong()
        val up = Regex("UP.*?recv_elapsed=(\\d+) ns").find(meta)!!.groupValues[1].toLong()
        return Fixture(readCsv("/fixtures/$name/acc.csv"), readCsv("/fixtures/$name/gyr.csv"), down, up)
    }

    /** 4.4 s 的慢滑：松手时手机基本已停（|v_end| 落在阈值内）→ 走裸积分，不顺延。 */
    @Test
    fun slowSlow_4p4s_matchesPython() {
        val f = load("20261011_090006")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)

        assertTrue("引擎应成功：${r.failReason}", r.ok)
        assertEquals(1625L, f.acc.size.toLong())
        assertEquals("DOWN 下标", 1091L, r.downIndex.toLong())
        assertEquals("UP 下标", 1553L, r.upIndex.toLong())
        assertEquals("终点（不该顺延）", 1553L, r.endIndex.toLong())
        assertFalse("不该顺延", r.extended)
        assertFalse("不该去趋势", r.detrended)

        // Python: slide_pipeline.py → ★ 距离 = 15.91 cm
        assertEquals("距离", 0.1590860109172647, r.distanceM, 1e-5)

        // 自校准的本机参数（这三项决定所有阈值，最容易在移植时写错）
        assertEquals("实测 |g|", 9.912948065697552, r.gravityMag, 1e-6)
        assertEquals("acc 噪声地板", 0.00818369183437402, r.accSigma, 1e-6)
        assertEquals("残余零偏 |b|", 0.0028311584305914764, r.biasBound, 1e-6)
    }

    /** 1.8 s 的快滑：松手时手机还有 2.0 cm/s → 触发顺延，终点从 #608 挪到 #632。 */
    @Test
    fun fastSlide_1p8s_extendsEndpointToRealStop() {
        val f = load("20261011_091811")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.upNs)

        assertTrue("引擎应成功：${r.failReason}", r.ok)
        assertEquals("DOWN 下标", 418L, r.downIndex.toLong())
        assertEquals("UP 下标（松手点）", 608L, r.upIndex.toLong())
        assertEquals("终点（顺延后）", 632L, r.endIndex.toLong())
        assertTrue("应发生顺延", r.extended)
        assertEquals("顺延 228 ms", 228.0, r.extendedMs, 1.0)
        assertFalse("顺延后 v_end 回到 0，不必再去趋势", r.detrended)

        // Python: 松手点只给 13.72 cm，顺延后 15.63 cm（真值 15 cm）
        assertEquals("距离", 0.15631377843131994, r.distanceM, 1e-5)
        assertEquals("顺延后的终点速度", -0.00022862142938964648, r.vEndFwd, 1e-6)
    }

    /** 窗口退化（按下即松手）必须被拒绝，而不是给一个看似合理的数字。 */
    @Test
    fun degenerateWindow_isRejected() {
        val f = load("20261011_091811")
        val r = SlideMeasure.measure(f.acc, f.gyr, f.downNs, f.downNs)
        assertFalse(r.ok)
        assertTrue("应给出原因", r.failReason != null)
    }

    /** 样本数不足必须被拒绝，而不是越界崩溃。 */
    @Test
    fun tooFewSamples_isRejected() {
        val f = load("20261011_091811")
        val r = SlideMeasure.measure(f.acc.take(10), f.gyr.take(10), f.downNs, f.upNs)
        assertFalse(r.ok)
    }
}
