# 把测距链路搬到设备端

> 之前 App 只采集，距离要 `adb pull` 回 PC 用 `sim/slide_pipeline.py` 算。
> 那就不是一把尺子，是一个数据记录仪。本文记录这一环怎么补的、以及怎么证明它搬对了。

---

## 1. 为什么必须搬

用户按一下按住式按钮、松手、1.5 秒后**屏幕上直接出数字** —— 这才是产品。
链路里的每一步都已经被真机验证过（见 `docs/boundary-design.md`、`docs/calibration.md`），
所以这一环的正确性判据只有一个：**设备端算出来的数，必须和 Python 算出来的数一样**。

---

## 2. 搬了什么

`app/src/main/java/com/thripleq/slideruler/probe/SlideMeasure.kt`，纯逻辑、零 Android 依赖，
与 `sim/pipeline.py` + `sim/slide_pipeline.py` 一一对应：

| 步骤 | Python | Kotlin |
|---|---|---|
| 时间基准 | `median(acc_ts 间隔)/1e9` | `median(dts)/1e9` |
| 边界对齐 | `locate(acc_ts, recv_elapsed_ns)` | `locate(acc, ns, n)` |
| 自校准 | `DeviceProfile`（DOWN 前 2 s 静止段） | `DeviceProfile`（同） |
| 姿态 | `solve_attitude`（互补滤波 α=0.05） | `qDelta` / `rotateInv` / `cross` |
| 去重力 | `remove_gravity(g_mag=prof.gravity_mag)` | `fw.z - profile.gravityMag` |
| 积分 | `integrate_vec`（三维梯形） | `integrateVec` |
| 距离 | `integrate_cond_zupt3d(k=5)` | `condZupt3d` |
| 顺延 | `stop_endpoint` | `stopEndpoint` |

**刻意没搬**的两项：

- **单姿态标度修正 `g_scale`**：`k̂ = ‖g‖/9.81` 在单姿态下混着零偏在重力方向的投影 `b·ĝ`，
  分不开（`b_z` 只要 0.05 m/s² 就装成 +0.5% 的标度）。
  参数保留只为逐位对齐，**默认 1.0 = 不修正**。详见 `pipeline.Params.use_scale_corr`。
- **静止检测器自动定位窗口**：按钮边界已经够准，检测器只当对照用，设备端不需要。

---

## 3. 怎么证明搬对了

`app/src/test/java/.../SlideMeasureTest.kt`，夹具是两条**真实采集**（连同 `meta.txt` 里的
触摸时刻），期望值由 `sim/slide_pipeline.py` 在同一份 CSV 上产出。

```
cd android && ./gradlew testDebugUnitTest
```

断言的不只是距离，还有**下标和分支选择** —— 那才是判据本身：两套实现都可能跑出接近
15 cm 的数字，但"用不用去趋势 / 要不要顺延终点"差一点点，换一条采集就会崩。

| 夹具 | DOWN→UP | 终点 | 分支 | 距离 | Python |
|---|---|---|---|---|---|
| `20261011_090006`（4.4 s 慢滑） | #1091→#1553 | #1553 | 裸积分、不顺延 | **15.909 cm** | 15.909 cm |
| `20261011_091811`（1.8 s 快滑） | #418→#608 | **#632** | 松手时还在动 → 顺延 228 ms | **15.631 cm** | 15.631 cm |

容差 `1e-5 m`（0.001 cm）—— 比界面显示的 0.01 cm 精度还紧 10 倍，只留给浮点求和顺序的差异
（Python 的 `fmean`/`pstdev` 用精确累加，Kotlin 用朴素求和）。

另外两条负例：`DOWN == UP` 的退化窗口、样本数不足，都必须**被拒绝**而不是给一个看似合理的数字。

---

## 4. 接入界面时唯一要小心的地方

`stop_endpoint` 的搜索上限是 **1.2 s**：松手时手机还在动的话，终点是算法在松手**之后**
才找到的，缓冲区里得先有那段数据。所以 `MainActivity.SLIDE_TAIL_MS` 必须 ≥ 1.2 s
（现取 1500 ms，留 0.3 s 余量）。

换算成用户看到的一条规则：**松手后保持手机不动 1.5 s 才出数**。提前把它拿起来，
停稳区找不到，算法就退回"松手点"当终点 —— 那会少算 `½·v_end·T`（本次约 1.9 cm）。

这一条也顺带修正了旧引导里的一句错话。以前写的是「灯是红的就松手 = 本次作废」，
现在**不成立了**：红着松手也能救回来，只要松手后手机在 1.5 s 内真的停下。
