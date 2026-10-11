// Top-level build file. AGP 9 provides built-in Kotlin, so `kotlin-android` is
// not applied; only the Compose compiler plugin is applied per-module.
// 版本组合与同机的 cirro 工程一致（已验证可编译），避免踩 AGP/Kotlin 版本坑。
plugins {
    alias(libs.plugins.android.application) apply false
    alias(libs.plugins.compose) apply false
}
