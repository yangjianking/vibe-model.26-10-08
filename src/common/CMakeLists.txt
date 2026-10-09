# common 模块：公共基础设施（类型、执行策略、配置、数学工具、性能基准）
# 大部分为 header-only，仅 config.cpp 与 performance.cpp 需要编译。

add_library(cs_common_impl OBJECT
    src/common/config.cpp
    src/common/performance.cpp
)
target_link_libraries(cs_common_impl PUBLIC cs_common yaml-cpp)

# 汇总为公共库别名（供其他模块链接）
add_library(cubed_sph::common ALIAS cs_common_impl)
