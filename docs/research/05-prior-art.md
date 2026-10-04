# 05 · 对口先例与竞品

## 最接近的先例：Wojtek120/IMU-velocity-and-displacement-measurements

**github.com/Wojtek120/IMU-velocity-and-displacement-measurements**（MIT，108★，README 已抓 `_raw/p3_wojtek.json`）

硕士论文项目：举重训练参数测量（杠铃的速 度/位移），架构与我们 spec **逐条对上**：

| 他的模块 | 我们 spec 的对应 |
|----------|------------------|
| MPU9250 + Arduino + 蓝牙 | 手机内置 IMU |
| Madgwick 滤波 | 姿态解算（我们候选之一） |
| Compensation of gravitational acceleration | 去重力 |
| **Zero Velocity Update (ZVU)** | 零速修正（连名字都一样） |
| 积分得速度 → 积分得位移 | 梯形二次积分 |
| 阈值 + 样本数窗的静止检测，参数"实验选定" | 我们的静止检测设计 |
| 已知起止静止 → 速度漂移线性扣除 | 线性漂移校正 |

关键实证：

- ZVU 静止检测 = "加速度与角速度不超过阈值并持续给定样本数"（阈值与窗长实验选定，无公式化——**这暗示我们要靠真机标定，文献没现成数**）
- 测试方法：机械臂带动给定距离/速度。**平均误差 ~12%**，作者归因于 $2 级便宜传感器
- 手机侧 Android App（蓝牙收数 + SQLite 存历史）——历史记录功能的先例

→ 这是最强的可行性锚点：**便宜的 MEMS + 同样算法管线 = 12% 误差**，正好落在 spec 预期带（5~15%）。我们条件更好（手机 IMU 好于 $2 模组、滑动比举重干净），有理由瞄准带内偏优。

## 教育级经典：CH Robotics

*Accelerometer, Position, and Orientation*（chrobotics.com/library/accel-position-velocity）——把"为什么双重积分短距离就崩"讲得最清楚的一篇科普，Arduino 论坛的经典引用。**开工前值得全读**（未抓取正文，短文，直接网页读）。

## 失败共识：Arduino MPU6050 论坛

forum.arduino.cc/t/237391（2014 起，抓取件 `_raw/p4_arduino.json`）——爱好者反复验证"便宜 IMU 双积分 hopeless"，被引 CH Robotics 收尾。价值：所有"我拿 MPU6050 积分测距"的野路子都会撞同一堵墙，我们靠的是 ZUPT + 短时窗 + 一维约束翻墙，不是硬积分。

## 车辆速度估计（同思路异场景）

Ustun, *Speed Estimation Using Smartphone Accelerometer Data*（Old Dominion Uni，PDF 直链）——用**起止点校准**（已知初速/末速为零）限制积分误差，与我们 ZUPT 首尾归零同构。参考其误差分析框架。

## iOS / 同形态竞品扫描

- CMPedometer 步行距离：公开数据 ~87.5% 准确率、±0.5 m 级（Watkins app 测试，ResearchGate）——步数×步长法，不是积分法
- **AR 测距（Measure app / 各家 AR ruler）**：摄像头方案，正是 spec 里要绕开的（弱光/纯色失效）
- **没有任何"贴边滑动测距"形态的产品或论文**——形态本身是空白，差异化成立
- GitHub 上 IMU odometry 项目多为机器人/VIO 辅助（`_raw/q6_projects.json`），无一维测距形态

## 综合定位

技术上：等价问题（IMU 短时位移 + 零速边界）已被充分探索，算法不必发明，**照抄成熟的**；
产品上：交互形态（贴边滑一滑，测个桌子长度）没人做——项目的价值全在"把成熟算法放进一个新交互形态，并把它调到可用精度"。
