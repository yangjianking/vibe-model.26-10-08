# =============================================================================
#  Cubed-Sphere NWP Model  — 顶层构建脚本
# =============================================================================
#  项目：基于立方球网格（Cubed-Sphere Grid）的全球非静力数值天气预报模式框架。
#
#  设计约束（遵循《NWP 模式开发 AI-Agent 提示词集》P0 铁律）：
#    - 密集计算内核使用 C++17；Python 仅用于后处理/配置/驱动。
#    - 并行模型：MPI（进程间）+ OpenMP（节点内）；GPU 经 Kokkos 抽象层接入。
#    - 精度三级别名 Real / ComputeReal / StateReal，CMake 选项切换。
#    - 所有可调参数经 YAML 注入，禁止硬编码。
# =============================================================================

cmake_minimum_required(VERSION 3.20)

# ----------------------------------------------------------------------------
# 项目元信息
# ----------------------------------------------------------------------------
project(CubedSphereNWP
        VERSION 0.1.0
        DESCRIPTION "Cubed-sphere non-hydrostatic atmospheric model framework"
        LANGUAGES CXX)

# 默认 C++17（可被上层覆盖为更高标准）
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 导出 compile_commands.json，供 IDE/静态分析使用
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ----------------------------------------------------------------------------
# 构建选项（精度 / 后端 / 模块开关）
# ----------------------------------------------------------------------------
option(CS_USE_KOKKOS         "通过 Kokkos 抽象层接入 GPU/众核后端" ON)
option(CS_ENABLE_MPI         "启用 MPI 分布式并行"                        ON)
option(CS_ENABLE_OPENMP      "启用 OpenMP 节点内并行"                    ON)
option(CS_ENABLE_NETCDF      "启用 NetCDF-4/HDF5 并行 I/O"               ON)
option(CS_ENABLE_TESTS       "构建单元/组件测试（Catch2）"                ON)
option(CS_ENABLE_DA          "预留 4D-Var 同化模块接口（本阶段占位）"      ON)
option(CS_BUILD_DOCS         "构建 Doxygen 文档"                          OFF)

# 精度模式：fp64（默认）| mixed（主体 fp32，几何/归约 fp64）
set(CS_PRECISION "fp64" CACHE STRING "计算精度模式：fp64 或 mixed")
set_property(CACHE CS_PRECISION PROPERTY STRINGS "fp64" "mixed")

# 垂直坐标：terrain-following（地形跟随高度）| pure-height（纯高度）
set(CS_VERTICAL "terrain-following" CACHE STRING "垂直坐标方案")
set_property(CACHE CS_VERTICAL PROPERTY STRINGS "terrain-following" "pure-height")

# ----------------------------------------------------------------------------
# 编译选项
# ----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "构建类型" FORCE)
endif()

# 零容忍警告（提示词 P9 硬性要求）
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(-w2)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
    add_compile_options(-Minform=inform)
endif()

include(cmake/CompilerOptions.cmake)

# ----------------------------------------------------------------------------
# 外部依赖查找
# ----------------------------------------------------------------------------
# yaml-cpp（配置解析，必需）
find_package(yaml-cpp QUIET)
if(NOT yaml-cpp_FOUND)
    # 退化到 FetchContent 自动拉取（pinned 版本见 cmake/Dependencies.cmake）
    include(cmake/Dependencies.cmake)
endif()

# MPI（可选）
if(CS_ENABLE_MPI)
    find_package(MPI COMPONENTS CXX)
    if(MPI_CXX_FOUND)
        message(STATUS "MPI 后端：${MPI_CXX_COMPILER}")
    endif()
endif()

# OpenMP（可选）
if(CS_ENABLE_OPENMP)
    find_package(OpenMP)
endif()

# NetCDF-4（可选，P4 并行 I/O 主格式）
if(CS_ENABLE_NETCDF)
    include(cmake/FindNetCDF.cmake)
endif()

# Kokkos（可选，P3 GPU/众核抽象）
if(CS_USE_KOKKOS)
    find_package(Kokkos QUIET)
    if(Kokkos_FOUND)
        message(STATUS "Kokkos 后端：已发现")
    else()
        message(STATUS "Kokkos 未发现，退化到自研模板执行策略层（同语义）")
    endif()
endif()

# ----------------------------------------------------------------------------
# 公共头文件目录
# ----------------------------------------------------------------------------
add_library(cs_common INTERFACE)
target_include_directories(cs_common INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(cs_common INTERFACE cxx_std_17)

if(CS_ENABLE_MPI AND MPI_CXX_FOUND)
    target_compile_definitions(cs_common INTERFACE CS_HAVE_MPI)
    target_include_directories(cs_common INTERFACE ${MPI_CXX_INCLUDE_DIRS})
    target_link_libraries(cs_common INTERFACE MPI::MPI_CXX)
endif()

if(CS_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(cs_common INTERFACE OpenMP::OpenMP_CXX)
endif()

# 精度别名宏定义（Real / ComputeReal / StateReal）
if(CS_PRECISION STREQUAL "mixed")
    target_compile_definitions(cs_common INTERFACE CS_PRECISION_MIXED)
endif()

# ----------------------------------------------------------------------------
# 子目录（各模块）
# ----------------------------------------------------------------------------
add_subdirectory(src/common)
add_subdirectory(src/parallel)
add_subdirectory(src/grid)
add_subdirectory(src/dynamics)
add_subdirectory(src/timeint)
add_subdirectory(src/physics)
add_subdirectory(src/io)
add_subdirectory(src/da)
add_subdirectory(src/diagnostics)
add_subdirectory(src/testcases)
add_subdirectory(src/driver)

# 测试
if(CS_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
