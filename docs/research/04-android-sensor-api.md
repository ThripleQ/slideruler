# 04 · Android 传感器 API 笔记

来源：developer.android.com 官方文档（*Motion sensors*，正文已抓 `_raw/android-motion-sensors.json`；*Reporting modes* from source.android.com，`_raw/android-reporting-modes.json`；*SensorEvent* API 参考，`_raw/sensorevent.json`）。抓取日 2026-10-05。

## 坐标系（官方原文实测规则）

设备平放桌面、自然朝向时：

- 左侧被推（设备向右移）→ **x 加速度为正**
- 底部被推（远离你）→ **y 为正**
- 向天空推 A m/s² → z = **A + 9.81**（静止时 z = +9.81）
- 陀螺：逆时针旋转（从正轴看向原点）为正；坐标系与加速度计一致

## 软件传感器陷阱（重要）

官方分类：`gravity`、`linear acceleration`、`rotation vector`、`significant motion`、`step counter` 等**是硬件或软件实现的派生传感器**，可用性和数据来源因设备而异（有的从加速度计+磁力计派生，有的从陀螺派生）。

→ **决策依据：`TYPE_LINEAR_ACCELERATION` 的去重力质量是厂商黑盒**（典型实现是低通估计重力，延迟与质量不可控）。spec 里"自己姿态投影去重力"的路线有官方文档背书。文档同时提醒：用原始加速度计时"you might have to implement low-pass and high-pass filters to eliminate gravitational forces"——正是我们要做的。

## 陀螺 rate limit（Android 12+）

targetSdk ≥ 31（API 31）时陀螺仪**受系统限速**。100 Hz 远低于限速线，不受影响；但如果将来想上 200+ Hz，需要 `HIGH_SAMPLING_RATE_SENSORS` 权限。→ manifest 里先留注释，不申请。

## 报告模式与采样

- 模式：continuous / on-change / one-shot / special（source.android.com *Reporting modes*）
- 我们要的加速度计/陀螺是 continuous 模式；用 `registerListener(listener, sensor, samplingPeriodUs)` 以**微秒**指定周期（100 Hz = 10000 µs）
- 官方最佳实践：指定能接受的最大 delay，设了就别中途改
- **`event.timestamp` 是纳秒、单调、适合做积分 dt**——用相邻事件时间差而非假设均匀间隔（论坛/文献都抱怨过采样不均匀，jitter 实际存在；这是选梯形积分 + 实测 dt 的原因）
- 传感器 FIFO batch：短时高率采样时留意 batching 会一次性吐历史数据（延迟观感差），必要时对 samplingPeriodUs 传 0 之外的显式值并禁用 batch FIFO（`SensorManager` 无直接 API，靠 rate hint；真机阶段验证）
- 功耗：官方说加速度计功耗约为其他运动传感器的 **1/10**——只在一维且姿态稳定时可考虑关陀螺省电，但第一版不做这种优化

## 采集实现备忘（为开工留）

- 后台测量若要保采样，前台 Service + partial wake（第一阶段屏亮测量，不需要）
- `Sensor.getResolution()` / `getMaximumRange()` 可做"传感器体检"展示（对应 03 号笔记的 Allan 方差体检思路）
- 多传感器时间对齐：accel 和 gyro 事件流独立回调，按各自 timestamp 对齐（最近邻或内插）——Niu 论文的课程设计也是这么处理的
