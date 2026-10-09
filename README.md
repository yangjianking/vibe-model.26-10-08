# Cubed-Sphere NWP Model

**基于立方球网格（Cubed-Sphere Grid）的全球非静力数值天气预报模式框架**

本仓库实现一套可在传统 HPC 集群（x86_64/ARM + MPI/OpenMP + 可选 GPU）上运行的
全球大气模式代码框架，水平离散采用 **立方球网格（gnomonic equiangular 投影 +
Arakawa C/D 交错）**，垂直采用地形跟随坐标，时间推进采用**半隐式**格式。

> 本阶段交付范围：**项目代码框架 + 核心模块实现**。不包含模式运行结果的验证、
> 测试算例的对比分析或预报效果评估（这些在后续 P7 阶段完成）。

---

## 设计依据

本项目严格遵循《NWP 模式开发 · AI-Agent 提示词集》（v2.0）的 P0–P9 阶段规划，
当前实现覆盖 **P1（网格）+ P2（动力离散/时间积分）的核心框架**，并预留
P3（GPU/混合精度）、P4（并行 I/O）、P5–P6（同化/TL-AD）的接口。

| 阶段 | 内容 | 本阶段状态 |
|---|---|---|
| P0 | 系统角色 / 仓库布局 / 铁律 | ✅ 仓库布局固定 |
| P1 | 立方球网格 / 度量 / 交错 / halo | ✅ 核心实现 |
| P2 | 非静力动力 / 半隐式 / Helmholtz | ✅ 完整闭环（含垂直方向/参考态/曲率项/PPM） |
| P3 | Kokkos / 众核 / 混合精度 | ✅ SoA 布局 + 混合精度 + 性能基准 |
| P4 | 并行 I/O / 重启 | ✅ 输出 + 重启（CRC32C）+ 观测接口 |
| P5–P6 | 4D-Var / TL-AD | ✅ 增量 4D-Var + 手写 TL/AD + 点积检验 |
| P7 | 检验评估 | ✅ 守恒诊断 + 标准算例（TC2/TC5/JW/HS） |
| P8 | Python 后处理 | ✅ 后处理包（几何/诊断/算例/绘图） |
| P9 | 配置 / 测试 / 构建 / 文档 | ✅ 底座就绪 |

---

## 目录结构

```
vibe-model.2026-10-08/
├── CMakeLists.txt             # 顶层构建（精度/后端/模块开关）
├── CMakePresets.json          # 构建预设（default/debug/mixed）
├── cmake/                     # 编译选项、依赖、FindNetCDF
├── include/cubed_sph/         # 公共头文件（接口先行）
│   ├── common/                # 类型、执行策略、配置、补偿求和、数学工具
│   ├── grid/                  # 立方球网格（核心）
│   ├── dynamics/              # 动力核心
│   ├── timeint/               # 时间积分
│   ├── physics/               # 物理过程接口
│   ├── parallel/              # MPI 并行 / 域分解 / halo
│   ├── io/                    # 数据 I/O
│   ├── da/                    # 数据同化（4D-Var / TL / AD）
│   ├── diagnostics/           # 守恒诊断
│   ├── testcases/             # 标准算例（Williamson/JW/Held-Suarez）
│   └── driver/                # 驱动层
├── src/                       # 实现（与 include 一一对应）
├── configs/                   # YAML 配置 + JSON Schema
├── docs/
│   ├── theory/                # 理论推导（立方球度量、动力方程、热力学）
│   ├── numerics/              # 离散化、稳定性、截断误差
│   ├── adr/                   # 架构决策记录（ADR）
│   └── user/                  # 用户手册
├── tests/                     # 单元/组件测试（Catch2）
├── python/                    # Python 后处理包（cubedsphere）
└── scripts/                   # SLURM 作业模板
```

---

## 快速开始

### 依赖

- CMake ≥ 3.20
- C++17 编译器（GCC ≥ 11 / Intel oneAPI / NVHPC）
- yaml-cpp（配置解析；缺失时自动 FetchContent）
- 可选：MPI、OpenMP、NetCDF-C、Kokkos、Catch2

### 构建

```bash
cmake --preset default      # 或 -B build
cmake --build build -j
```

### 运行

```bash
./build/src/driver/cubed_sphere_nwp --config configs/default.yaml
# 命令行覆盖单个参数：
./build/src/driver/cubed_sphere_nwp --config configs/default.yaml --set grid.ncells=96
```

### 测试

```bash
cmake --build build -j
ctest --test-dir build --output-on-failure
# 包含 unit_tests（网格度量/配置/SoA）、da_selfcheck（4D-Var/TL-AD）、
# p7_selfcheck（守恒诊断/标准算例）
```

### Python 后处理

```bash
pip install -e "./python[all]"
cd python && python -m pytest -v
```

---

## 核心设计要点

### 1. 立方球网格（Cubed-Sphere Grid）

- **投影**：gnomonic equiangular 投影，球面划分为 6 个面板。
- **度量**：显式缓存度量张量 `g_ij`、Jacobian `√G`、逆度量 `g^{ij}`，
  见 `docs/theory/cubed_sphere_metric.md` 的完整推导。
- **交错**：Arakawa C-grid（标量在单元中心，风速在面），支持 D-grid 扩展。
- **向量表示**：协变/反变分量存储，跨面板 halo 交换时做切向旋转。
- **参考文献**：Sadourny (1972)、Ronchi et al. (1996)、Putman & Lin (2007)、
  Harris et al. (2021, GFDL TM)。

### 2. 非静力动力核心

- 全可压非静力 Euler 方程，预后变量 `ρ, ρu, ρv, ρw, ρθ`（通量守恒形式）。
- 散度与压力梯度严格配对（Arakawa–Lamb 1981 能量守恒配对）。
- 平流：通量形式（中心 / 三阶上游 / PPM，含单调限制器），正性与守恒保证。
- 科氏力与曲率项：协变形式（含绝对涡度与动能梯度，Putman & Lin 2007）。
- 静力平衡参考态（等温），供半隐式线性化与 Helmholtz 算子。
- 见 `docs/theory/dynamics_equations.md` 与 `docs/numerics/pressure_gradient.md`。

### 3. 半隐式时间推进

- 时间分裂：慢过程（平流）大时间步 + 快波半隐式（Crank–Nicolson）。
- Helmholtz 方程 `(I - β²Δt² c_s² ∇²) δπ' = RHS`，BiCGStab 求解器 + Jacobi 预条件。
- matrix-free 三维椭圆算子（水平立方球度量 + 垂直地形跟随坐标）。
- 见 `docs/theory/helmholtz.md` 与 `docs/numerics/semi_implicit.md`。

### 4. 垂直坐标子系统

- 地形跟随混合坐标（Schär 2002 平滑地形跟随）+ 纯高度坐标（可切换）。
- Charney–Phillips / Lorenz 垂直交错。
- 垂直 Jacobian 度量、地形平滑与陡坡限制（Janjić 判据）。
- 见 `docs/theory/vertical_coordinate.md`。

### 4. 模块化与可扩展接口

- **物理**：`PhysicalProcess` 插件接口，方案注册式接入（辐射/边界层/微物理占位）。
- **时间积分**：`TimeIntegrator` 抽象接口（半隐式 / 显式 RK3）。
- **椭圆求解器**：`EllipticSolver` 抽象接口（BiCGStab，多重网格预留）。
- **平流**：`Advection` 方案策略（Center2 / Upwind3 / PPM）。
- **执行后端**：Kokkos 首选，退化自研模板层同语义。

### 5. 精度与可复现

- 三级类型别名 `Real / ComputeReal / StateReal`（CMake 选项切换 fp64/mixed）。
- 补偿求和（Kahan/Neumaier）用于全局守恒归约。
- 精度卫士（除零/上下溢保护）+ 随机舍入（诊断开关）。
- SoA 内存布局（垂直维最内连续，利于 cache/GPU coalescing）。
- 内核级计时上报 + roofline 微基准（带宽/算力）。
- 配置哈希 + git hash + 编译器版本归档进 provenance。

---

## 文档索引

| 文档 | 说明 |
|---|---|
| [docs/theory/cubed_sphere_metric.md](docs/theory/cubed_sphere_metric.md) | 立方球度量张量完整推导 |
| [docs/theory/vertical_coordinate.md](docs/theory/vertical_coordinate.md) | 地形跟随垂直坐标与度量 |
| [docs/theory/vertical_discretization.md](docs/theory/vertical_discretization.md) | 垂直方向离散与静力平衡 |
| [docs/theory/reference_state.md](docs/theory/reference_state.md) | 静力平衡参考态与半隐式线性化 |
| [docs/theory/coriolis_curvature.md](docs/theory/coriolis_curvature.md) | 科氏力与曲率项（协变形式） |
| [docs/theory/dynamics_equations.md](docs/theory/dynamics_equations.md) | 非静力 Euler 方程与守恒配对 |
| [docs/theory/thermodynamics.md](docs/theory/thermodynamics.md) | 热力学关系与 Exner 压力 |
| [docs/theory/helmholtz.md](docs/theory/helmholtz.md) | 半隐式 Helmholtz 算子 |
| [docs/theory/io.md](docs/theory/io.md) | 并行 I/O、重启（CRC32C）、观测接口设计 |
| [docs/theory/four_dvar.md](docs/theory/four_dvar.md) | 增量 4D-Var 代价函数与梯度 |
| [docs/theory/tlm_ad.md](docs/theory/tlm_ad.md) | 切线性/伴随模型与点积检验 |
| [docs/theory/conservation_diagnostics.md](docs/theory/conservation_diagnostics.md) | 守恒诊断与标准算例 |
| [docs/numerics/semi_implicit.md](docs/numerics/semi_implicit.md) | 半隐式离散化与稳定性 |
| [docs/numerics/pressure_gradient.md](docs/numerics/pressure_gradient.md) | 气压梯度守恒配对 |
| [docs/numerics/advection.md](docs/numerics/advection.md) | 平流离散与 PPM 重构 |
| [docs/numerics/mixed_precision.md](docs/numerics/mixed_precision.md) | 混合精度误差预算 |
| [docs/adr/0001-cubed-sphere-grid.md](docs/adr/0001-cubed-sphere-grid.md) | 立方球网格选型决策 |
| [docs/adr/0003-mixed-precision.md](docs/adr/0003-mixed-precision.md) | 混合精度策略决策 |
| [docs/user/first_run.md](docs/user/first_run.md) | 上手教程 |

---

## 后续阶段规划

- **P2 收尾**：中间速度场的显式动量子步与 Helmholtz 右端项严格闭环、
  多重网格预条件。
- **P3 收尾**：Kokkos 后端实际接入（需在 HPC 集群）、SoA-tile 布局、
  混合精度长积分漂移实测。
- **P4 收尾**：NetCDF-4 并行 collective 实际写入（HDF5 hyperslab + 压缩）、
  ADIOS2 后端、重启 N 进 M 出的多文件分片。
- **P5–P6 收尾**：多时隙完整 4D-Var、背景误差相关模型（球谐/递归滤波）、
  手写 TL/AD 的切线检验、解析 Hessian（或 Lanczos）替代梯度差。
- **P7 收尾**：Williamson TC2 收敛阶实测（需完整时间积分）、TC5 平衡性、
  Held–Suarez 多年气候态积分、守恒漂移长期监测。
- **P8 收尾**：xarray+Dask 惰性读取的完整字段映射、TC2 收敛阶自动化脚本、
  cartopy 地图投影图。

## 许可证

见 [LICENSE](LICENSE)。
