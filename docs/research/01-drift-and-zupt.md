# 01 · 双重积分漂移与 ZUPT

## 误差的本质（为什么纯积分必炸）

加速度计双重积分的误差传播规则（El-Sheimy 的经典结论，被广泛引用）：

- **零偏（bias）**：速度误差 **线性** 增长，位置误差 **平方** 增长
  - δa = 0.05 m/s² 的恒定零偏，3 秒后位置误差 = ½·δa·T² = **0.225 m**——已经吃满预算
- **白噪声**：速度误差 ~√T 增长，位置误差 ~T^1.5 增长（量级换算见 03 号笔记）
- 去重力不干净（姿态误差 θ 使重力分量泄漏）：
  δa = g·sin(θ) ≈ 9.81·θ；**1° 姿态误差 ≈ 0.17 m/s² 假加速度**——比传感器零偏大一个量级
  → 这就是为什么姿态解算精度是一维测距的命门，而不是传感器本身

社区共识（Arduino MPU6050 论坛，2014，原帖抓取件 `_raw/arduino-mpu6050-forum.json`）：

> "In principle it is possible to estimate the velocity and distance from the acceleration by double integration, but in practice with hobby-type accelerometers, this works for only very short distances before the estimate becomes hopelessly inaccurate." —— 指向 CH Robotics 的经典解释文。⚠️ 抓取时（2026-10-05）`chrobotics.com` 域名已易主（返回无关内容），原文不可达；其核心论证（双重积分误差平方爆炸、只有短时+边界约束可用）与本笔记推导一致，已转述于上

## ZUPT（零速修正）文献

| 文献 | 要点 |
|------|------|
| Foxlin 2005, *Pedestrian Tracking with Shoe-Mounted Inertial Sensors*, IEEE Computer Graphics & Applications | ZUPT 鼻祖级应用：鞋挂 IMU，每步触地瞬间强制速度归零，误差不随总时长累积。**关键洞察：ZUPT 把误差增长从"随总时间"变成"随步间时长"** |
| Zhang et al. 2022, *Deep neural network-based adaptive zero-velocity detector* (IET) | 用 DNN 做自适应静止检测（阈值法在非理想步态下失效时的改进路线） |
| *Understanding the performance of zero velocity updates in pedestrian navigation* (Nottingham) | 系统分析 ZUPT 何时有效、何时失效 |
| Li et al. 2012, *A robust pedestrian navigation algorithm with low cost IMU* (IPIN) | 低成本 IMU + ZUPT 的工程实现参考 |
| Hölzke et al. 2019, *Low-complexity online correction and calibration* (Taylor & Francis) | **与我们最像**：ZUPT 算单步步长 + Madgwick 滤波，低复杂度在线校正（正文 paywall，只有摘要级信息） |

## 对本项目的特殊有利条件

文献里 ZUPT 都用在**步态**场景——"触地瞬间速度为零"只有几十毫秒、且伴随剧烈冲击振动，检测难。
我们的场景是**贴桌直线滑动**：

- 起点和终点有**真正的长时间静止**（用户主动停 1~2 秒）→ 静止检测极可靠
- 没有冲击振动，信号干净
- 一维直线约束 + 0.2~2 m / 1~3 s 的短时窗

即：**我们是 ZUPT 的理想工况，比步态场景容易得多**。这是 spec 精度预期的理论根基。

## 待验证

- 静止检测阈值（加速度模长偏差 + 陀螺范数）的取值范围——文献没有现成数，仿真 + 真机标定
- 中途犹豫停顿（停 0.5s 又继续）会不会把一次测量切成两段——交互设计上要处理
