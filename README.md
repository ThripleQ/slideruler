# slideruler

**基于手机 IMU 的滑动式一维测距工具。**

用户把手机一角贴住桌沿，从一端平稳滑到另一端，App 只靠加速度计 + 陀螺仪（不调用摄像头、不用 AR、不连外部硬件）算出滑动距离——即被测物体的长度。

面向 0.2 ~ 2 m 的直线边缘测量（桌子、柜子、板材），目标误差 5% ~ 15%。

```
校准(静止2s) → 100Hz采集 → 姿态解算 → 去重力 → 一维约束 → ZUPT零速修正 → 二次积分 → 距离
```

## 当前状态：仿真验证完成 → 采集层 + 体检 A + 按住式边界落地

- **想弄懂这个项目的技术问题** → 从 [`docs/research/00-why-this-is-hard.md`](docs/research/00-why-this-is-hard.md) 开始（为什么难、靠什么救、精度算账）
- **想知道怎么一步步实现**（采什么数据、数据有什么缺陷、每步怎么验证）→ [`docs/implementation-guide.md`](docs/implementation-guide.md)（九阶段小白向全流程指南）
- 调研笔记总索引：[docs/research/README.md](docs/research/README.md)——从"双重积分为什么会漂"到"手机 IMU 实测噪声数字"到"最接近的先例（同管线实测 12% 误差）"
- **仿真实验记录** → [`sim/README.md`](sim/README.md)：用已知真值的合成数据跑完整管线，验证算法无结构性错误，并把精度瓶颈定位到「积分起点边界」；含指南 §7 的三处勘误。
- **真机采集 / 体检 A / 滑动测量** → [`android/README.md`](android/README.md)：Android 工程，把原始 IMU（加速度计 + 陀螺，100 Hz）落盘成 CSV，在手机上直接给这台设备做静态体检，并支持**按住式滑动测量**（触摸事件当积分边界）。
- **体检 A 实测报告** → [`docs/health-a.md`](docs/health-a.md)：K40（lsm6dso）的实测零偏/噪声/时间戳质量，以及"去重力基准吸收零偏"的实证。
- **算法校准设计** → [`docs/calibration.md`](docs/calibration.md)：为什么必须做「每会话自校准」（时间基单项代价 +11%）、单姿态标度修正为什么是陷阱、多姿态标定怎么做。
- **边界设计** → [`docs/boundary-design.md`](docs/boundary-design.md)：把"算法精确检测边界"换成"按住式手动边界 + 算法审计"。窗口边界只要落在静止段内，误差就与深度无关（0.12%，真机已验证）。

下一阶段：真机做一次真实滑动（验证端到端 + 手指按压干扰）→ 审计结果进 App UI → 姿态/去重力/PCA 进 App + 回放 UI。

```
slideruler/
├── docs/            资料包 + 实现指南 + 体检报告 + 校准设计 + 边界设计
├── sim/             合成数据仿真 + 真机数据分析（纯标准库 Python）
└── android/         Android 采集层 + 体检 A + 按住式滑动测量（Kotlin / Compose）
```
