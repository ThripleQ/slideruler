# slideruler 资料包 —— 基于 IMU 的滑动式一维测距

搜集时间：2026-10-04 ~ 10-05（web 搜索 + 页面抓取，原始件存 `_raw/`）。

## 结论速览

| # | 结论 | 对项目意味着什么 |
|---|------|------------------|
| 1 | 加速度双重积分的误差：零偏 → 速度线性增长 → 位置平方增长；白噪声同样不可忽略 | 纯积分必炸，误差控制全靠 ZUPT / 一维约束 / 短时窗（spec 的三板斧方向正确） |
| 2 | 同类先例（举重位移测量，$2 级 IMU + Madgwick + ZVU）机械臂实测误差 ~12% | spec 的 5%~15% 预期有实证支撑，不是空想 |
| 3 | 手机 IMU 实测：加速度白噪声 ~6.5 m/s/√h，零偏不稳定性 ~30 mGal | 仿真参数有了实数来源（见 03 号笔记的量级换算表） |
| 4 | Android 的 `TYPE_LINEAR_ACCELERATION` 是软件派生传感器，实现质量因厂商而异 | 去重力自己做（姿态投影），不吃系统黑盒——与 spec 一致 |
| 5 | 没找到"贴边滑动测长"这个具体形态的已有产品/论文 | 差异化成立；但等价问题（IMU 短时位移）文献充分，不必重新发明算法 |

## 目录

- [00-为什么难.md](00-why-this-is-hard.md) —— 技术问题入门：误差三杀手、ZUPT 为什么是解药、精度算账（**小白先读这份**）
- [01-双重积分漂移与ZUPT.md](01-drift-and-zupt.md) —— 误差本质 + 零速修正文献
- [02-姿态滤波选型.md](02-attitude-filters.md) —— 互补 / Madgwick / Mahony 对比与选型建议
- [03-手机IMU特性实测.md](03-smartphone-imu-characteristics.md) —— 噪声/零偏实测数字 + 对本项目时窗的量级换算
- [04-Android传感器API笔记.md](04-android-sensor-api.md) —— 官方文档要点：坐标系、软件传感器、报告模式
- [05-对口先例与竞品.md](05-prior-art.md) —— Wojtek 项目（最接近的先例）、CH Robotics、iOS 侧情况
- [`../implementation-guide.md`](../implementation-guide.md) —— **实现全流程指南**：九阶段每步做什么/数据长什么样/缺陷/怎么验证（小白主文档）
- `_raw/` —— 搜索结果与页面抓取的原始 JSON（可追溯）

## 资料缺口（如实记录）

1. Madgwick 2010 原始报告 PDF：x-io 官链 404，wayback 与 archive.org 镜像均未取得——**非阻塞**，算法细节有 arduino-libraries/MadgwickAHRS 官方源码可读（比论文更实用）
2. Hölzke 2019（ZUPT 步长 + Madgwick 低复杂度校正）正文在 paywall 后，只有摘要
3. CH Robotics *accel-position-velocity*：抓取时域名已易主（返回无关内容），原文不可达；核心论证已在 01 号笔记转述
4. 手机贴桌滑动姿态晃动的实测数据（决定一维 PCA 假设的鲁棒性）——**只能自己测**，属于真机阶段
