# android —— 采集层 + 体检 A + 滑动测量

> 对应 `docs/implementation-guide.md` 的实现顺序第 1~2 步：**先有数据工厂，再把真机传感器摸清**。
> 现阶段这个 App 做三件事：**把原始 IMU 落盘**、**给这台手机做静态体检**、
> 以及**按住式滑动测量**（用触摸事件当积分边界，见 `docs/boundary-design.md`）。
> 算法（姿态/去重力/PCA/ZUPT/积分）**不在这里**——按指南 §0 的架构决定，采集与计算分离，
> 算法在测量结束后对 CSV 离线跑（先在工作区的 `sim/` 里调，之后搬进 App）。

## 1. 编译与安装

环境（本机已验证）：Android Studio 自带 JBR（OpenJDK 25）、Android SDK `platforms/android-37.0`、
`build-tools/36.0.0`、Gradle 9.5.0（wrapper 指向的发行版已在本地缓存）。
版本组合（AGP 9.3.1 / Kotlin 2.4.10 / compileSdk 37 / minSdk 26）与同机已能编译的 `cirro` 工程一致。

```bash
cd slideruler/android
export JAVA_HOME="/c/Program Files/Android/Android Studio/jbr"   # 或直接在 Android Studio 里打开本目录
./gradlew :app:assembleDebug
./gradlew :app:installDebug          # 需要已连接设备（adb devices 能看到）
```

或直接 `adb install -r app/build/outputs/apk/debug/app-debug.apk`。

打开方式：Android Studio → Open → 选 `slideruler/android`（**不是**仓库根目录，根目录没有 Gradle 工程）。

## 2. 体检 A 怎么做

1. 手机**平放**在桌面上，拿走上面所有东西（壳、支架、绳子）。
2. 点「**体检 A · 静止 30 s**」，**不要碰桌子，也不要说话**，等倒计时走完（屏幕采集期不会灭）。
3. 自动结束 → 屏上直接给出报告，同时把 CSV 与档案写进手机。

> 这一步验的是 `docs/implementation-guide.md` §2 的四条：模长 / 零偏 / 噪声 / 间隔直方图。
> 报告里每一项后面带 ✓ △ ✗，判读规则就写在报告里，不用回来翻文档。

## 3. 滑动测量（按住式边界）怎么做

点「**滑动测量（按住式边界）**」进入。它会**立刻开始连续采集**，等你按住大按钮来划定边界。

1. 先把手机放好，**静置 2 秒以上**（这段用来做自校准：陀螺零偏、重力参考、噪声地板）。
2. **按住**大按钮 → 手指别抬 → **停顿约 0.5 秒**（这时候手机还是静止的）。
3. 推动手机到目标位置（手和手机一起走，手指始终压在屏幕上别松）。
4. **推到位后停稳**，再**松手**。
5. 松手后再采 0.7 秒自动结束并落盘。

**为什么是"按住 → 停顿 → 推 → 停稳 → 松手"这五个动作**：边界只要落在**静止段**里，
误差就与"落得多深"无关（`docs/boundary-design.md` §2 的物理推导 + §3 的实测）。
所以"停顿"和"停稳"这两步是刻意留给安全余量的；真正怕的是**按下就推**或**没停就松手**
——那会把积分窗口切进运动段，误差直接从 0.1% 跳到 −15%。

> 采集后 App 会做一次**审计**：检查 DOWN/UP 时刻附近 ±0.3 s 的线性加速度是否在静止阈值内。
> 审计不通过就提示重测。这是把"算法必须精确定位边界"换成"用户必须按对 + 算法验真"。

用 adb 模拟一次（不开真机也能验证机制）：`adb shell input swipe X Y X Y 3000`
（同点 swipe 就是长按 3 秒）。

## 4. 报告怎么读

| 项 | 合格线 | 说明 |
|---|---|---|
| ‖a‖ 均值 | 偏离 9.81 < 0.10 | 偏差 > 0.30 说明这颗加速度计有问题，先别往下做 |
| 三轴零偏 | ——（记录下来） | 这是**这台机器**的零偏，进「机型档案」，别用文档里的别人的数 |
| 三轴噪声 σ | ≤ 0.05 正常，> 0.08 难做 | 指南说 0.01~0.05 都正常；这条决定 `sim/` 里 `acc_noise` 该改成多少 |
| 实测采样率 | 接近 100 Hz | 有的机器你请求 100 Hz 它给 98 Hz，正常 |
| 间隔 P95 / max | P95 ≤ 13 ms | 尾巴长 = 调度不稳；**无论怎样都要按实测 dt 积分，永远不要写死 0.01** |
| 丢样 (dt>50 ms) | 0 次 | 非 0 说明有断流，算法层必须把这一段标记出来 |

## 5. 数据在哪、怎么取回

写在 app 专属外部目录（不需要任何存储权限）：

```
/sdcard/Android/data/com.thripleq.slideruler.probe/files/slideruler/<yyyyMMdd_HHmmss>/
    acc.csv      timestamp_ns,ax,ay,az      加速度计（m/s²，比力）
    gyr.csv      timestamp_ns,gx,gy,gz      陀螺（rad/s）
    meta.txt     设备/传感器档案 + 体检报告 + 触摸事件
```

滑动测量的 `meta.txt` 里会多一段触摸事件（**这是按钮边界的原始记录**）：

```
触摸事件（按住式边界）
  DOWN event=<uptimeMillis> ms  recv_uptime=<ms> ms  recv_elapsed=<ns> ns  相对首样本=<ms> ms  输入延迟=<ms> ms
  UP   event=...                recv_uptime=...      recv_elapsed=...        相对首样本=...        输入延迟=... 
```

三个时间戳的分工：`recv_elapsed_ns` **与 `acc.csv` 同基准**，可直接对齐；
`event` 与 `recv_uptime` 同为 uptime 基准，两者之差就是**输入延迟**。

取回（USB 连着时）：

```bash
ADB="$LOCALAPPDATA/Android/Sdk/platform-tools/adb.exe"
"$ADB" pull /sdcard/Android/data/com.thripleq.slideruler.probe/files/slideruler/ ./
```

两条流**分开成两个文件**是刻意的：它们的回调时刻本来就不重合（加速度计 12:00:00.010，
陀螺 12:00:00.013），时间对齐是算法层的事，别在采集层假装它们对齐了。

## 6. 采集层的几个刻意决定

- **只用两个原始硬件传感器**（`TYPE_ACCELEROMETER` / `TYPE_GYROSCOPE`）。
  不用 `TYPE_LINEAR_ACCELERATION` / `TYPE_GRAVITY`——厂商黑盒派生量，质量不可控（指南 §2 第 6 条）。
- **请求 100 Hz**（`registerListener(..., 10_000)` 微秒）。想测这台机器的上限，把
  `ImuRecorder.DEFAULT_PERIOD_US` 改成 `SensorManager.SENSOR_DELAY_FASTEST` 再跑一次即可
  （已声明 `HIGH_SAMPLING_RATE_SENSORS`，API 31+ 请求 >200 Hz 需要它）。
- **时间戳只用 `event.timestamp`**（纳秒、单调、同设备所有传感器共用一个时钟，官方保证）。
- 缓冲上限 30 万样本（100 Hz 下约 50 分钟），纯粹防 OOM。

## 7. 现在还没做的（下一步）

- **回放 UI**（指南 §10）：选历史 CSV → 在手机上跑完整管线 → 出体检 B~E。目前管线在工作区 `sim/` 里跑，
  真机版入口是 `sim/slide_pipeline.py`。
- **审计结果进 UI**：现在审计只在离线脚本里算，还没搬到 App 的结果卡片上。
- CSV 的**分享导出**（FileProvider）——目前靠 `adb pull` 够用。
- 采集期的**实时曲线**（指南说的"另一条轻量路径"）。

## 8. 与 sim/ 的关系

`sim/README.md` 里那份敏感性扫描的所有传感器参数（`acc_noise`、`acc_bias`、`gyro_bias`、
`gyro_instab`）原本是**合成数**。体检 A 跑完（`docs/health-a.md`）已经换成这台
**Xiaomi M2012K11AC（Redmi K40）** 的实测值，`sim/` 的默认参数按它更新过。
换机时要重做一遍体检 A，并按新机数字更新 `sim/` 的 `SensorCfg`。
