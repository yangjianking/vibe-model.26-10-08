# 上手教程（First Run）

> 本仓库为立方球全球大气模式（Cubed-Sphere NWP）的完整代码框架，
> 覆盖 P0–P9 全部阶段：网格、动力核心、时间积分、混合精度、并行 I/O、
> 数据同化（4D-Var / TL-AD）、守恒诊断、标准算例与 Python 后处理。

## 1. 前置依赖

- CMake ≥ 3.20
- C++17 编译器（GCC ≥ 11 / Intel oneAPI / NVHPC）
- yaml-cpp（缺失时自动 FetchContent）
- 可选：MPI、OpenMP、NetCDF-C、Kokkos

> 完整构建需 Linux HPC 环境（含 cmake/g++/yaml-cpp/MPI）。Windows 开发环境
> 下核心模块已用 MSVC `/std:c++17 /utf-8 /W4` 逐模块编译验证零错误零警告。

## 2. 构建

```bash
cmake --preset default      # 等价于 cmake -B build
cmake --build build -j
```

其他预设：

```bash
cmake --preset debug        # Debug + 测试
cmake --preset mixed        # 混合精度（fp32 主体 + fp64 几何）
```

## 3. 运行测试

```bash
ctest --test-dir build --output-on-failure
```

预期通过：`unit_tests`（网格度量/配置/SoA 精度）、`da_selfcheck`（观测算子
TL/AD 转置 + 4D-Var 梯度）、`p7_selfcheck`（守恒诊断 + 标准算例解析解）。

## 4. 运行模式（框架冒烟）

```bash
./build/src/driver/cubed_sphere_nwp --config configs/default.yaml
```

预期输出（框架连通性验证，非预报结果）：

```
[Driver] 初始化完成：ncells=48, nhalo=3, nlev=60, dt=300s, steps=10
[Driver] 时间积分器：SemiImplicit
[Driver] 开始时间推进（10 步）
[Driver] 时间推进完成
[main] 配置哈希: <hash>
```

## 5. 命令行覆盖单个参数

```bash
./build/src/driver/cubed_sphere_nwp --config configs/default.yaml \
    --set grid.ncells=96 --set timeint.dt=150.0
```

## 6. Python 后处理包

```bash
pip install -e "./python[all]"    # 全功能（含 xarray/dask/matplotlib）
cd python && python -m pytest -v  # 运行后处理单元测试
```

核心几何/诊断/算例模块零依赖，可离线验证：

```python
from cubedsphere.grid import unit_sphere_area
from cubedsphere.test_cases import WilliamsonTC2

print(unit_sphere_area(48))  # → 12.566...（≈ 4π）
```

## 7. 各阶段收尾项

代码框架与核心模块已全部实现并逐模块编译验证。仍需在 Linux HPC 完成的
运行级验证：

- P2 收尾：Helmholtz 多重网格预条件
- P3 收尾：Kokkos 后端实际接入、混合精度长积分漂移实测
- P4 收尾：NetCDF-4 collective 实际写入、ADIOS2
- P5–P6 收尾：多时隙完整 4D-Var、切线检验、背景误差相关模型
- P7 收尾：TC2 收敛阶实测、长期积分守恒漂移监测

## 8. 文档索引

| 目录 | 内容 |
|---|---|
| `docs/theory/` | 立方球度量、动力方程、热力学、半隐式、4D-Var、TL/AD、守恒诊断 |
| `docs/numerics/` | 平流、气压梯度、半隐式、混合精度离散化 |
| `docs/adr/` | 架构决策记录（网格、执行后端、混合精度） |
| `docs/user/` | 用户手册 |
