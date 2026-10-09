# =============================================================================
#  第三方依赖 pinned 版本（superbuild）
# =============================================================================
# 说明：网络不可用或未安装 yaml-cpp 时，通过 FetchContent 拉取固定版本，
#       保证可复现构建。生产 HPC 集群建议改用 Spack 环境（见 spack.yaml）。

include(FetchContent)

# yaml-cpp —— 配置解析（必需依赖）
set(YAML_CPP_TAG "0.8.0")
FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG        ${YAML_CPP_TAG}
    GIT_SHALLOW    TRUE
)
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(yaml-cpp)

# Catch2 —— 单元测试框架（仅测试目标需要）
if(CS_ENABLE_TESTS)
    set(CATCH2_TAG "v3.5.0")
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        ${CATCH2_TAG}
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(Catch2)
endif()

# 说明：NetCDF-C/C++、HDF5、Kokkos、METIS、ADIOS2 等重型依赖在 HPC 集群上
#       通常由系统模块或 Spack 提供，此处不默认拉取，由 FindXxx.cmake 探测。
