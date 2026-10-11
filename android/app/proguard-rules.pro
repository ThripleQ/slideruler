# 采集/体检工具没有反射、序列化、JNI，默认规则即可。
# 保留行号便于真机崩溃栈定位（见 android/README.md）。
-keepattributes SourceFile,LineNumberTable
-renamesourcefileattribute SourceFile
