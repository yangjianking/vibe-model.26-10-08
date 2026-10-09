# 未完成工作清单（TODO / Backlog）

> 本文档记录立方球网格 NWP 模式项目在 P0–P9 全部框架完成后的**剩余工作**。
> 代码框架（34 头文件 + 25 源文件 + 13 理论文档 + 7 Python 模块）已 100% 完成，
> 下列各项为「接口已就绪、实现待填充」的实质功能，或需 Linux HPC 环境的运行级验证。
>
> 最后更新：2026-10-08

---

## 一、必须依赖 Linux HPC 环境（当前 Windows 无法验证）

这些工作依赖 cmake/g++/yaml-cpp/MPI/NetCDF/Kokkos 等，当前 Windows 开发环境缺失，
需部署到 Linux HPC 集群完成。

| # | 工作项 | 位置 | 依赖 | 说明 |
|---|---|---|---|---|
| 1 | 完整 CMake 链接构建 | 全项目 | cmake/g++/yaml-cpp/MPI | `config.cpp` 依赖 yaml-cpp、`main.cpp` 需 MPI；当前仅逐模块编译到 obj，未做最终链接 |
| 2 | 跨进程 halo 交换 | `src/parallel/parallel.cpp` | MPI | `HaloExchange` 的 `MPI_Isend/Irecv` + 双缓冲、向量场交换 + 跨面板旋转（`panel_rotation`）均为 TODO |
| 3 | NetCDF-4 并行 collective 写入 | `src/io/io.cpp` | NetCDF/HDF5 | `nc_create` + 维度定义 + UGRID 元数据、`nc_put_vara_*` 并行 collective 均为 TODO，当前是二进制降级占位 |
| 4 | Kokkos 后端实际接入 | `common/execution_policy.hpp` | Kokkos/GPU | 当前用自研模板执行策略层退化路径，未接真实 Kokkos |

---

## 二、接口已就绪但实现为「占位/骨架」

这些功能接口完整、可编译、可链接，但核心算法是占位实现，需真实物理/数值工作填充。

| # | 工作项 | 位置 | 说明 |
|---|---|---|---|
| 5 | 物理参数化方案 | `src/physics/physics.cpp` | 辐射(RRTMG)、边界层(MYNN/YSU)、微物理(Thompson/WSM6) 全是 Dummy 占位，加热率恒为 0 |
| 6 | Helmholtz 多重网格预条件 | `src/timeint/` | 当前仅 BiCGStab + Jacobi 预条件，多重网格是 P2 收尾项 |
| 7 | 手写 TL/AD 切线检验 | `src/da/tlm_ad.cpp` | 已有 `FiniteDifferenceTLM` 参照，但未做切线检验（点积检验仅完成数学验证） |
| 8 | 背景误差相关模型 | `src/da/` | 背景误差协方差 B 目前为对角，球谐谱/递归滤波未实现 |
| 9 | 多时隙完整 4D-Var | `src/da/four_dvar.cpp` | 当前为单时隙接口联调骨架，多时隙 + 完整时间窗组织未做 |
| 10 | 跨面板 halo 交换接入平流 | `src/dynamics/advection.cpp` | 两阶段通量配对已守恒（实测机器精度），但全球完整守恒还需跨面板 `HaloExchange` 使面板边界通量一致 |

---

## 三、运行级验证（需完整时间积分跑通）

| # | 工作项 | 说明 |
|---|---|---|
| 11 | TC2 收敛阶全模式实测 | 已在 Python 侧完成平流算子收敛阶（Center2=2.0 / Upwind3=3.0 / PPM=2.35），未在完整模式（含时间积分器）下跑 |
| 12 | TC5 罗斯贝波平衡性 | 地形罗斯贝波长时间积分平衡性未验证 |
| 13 | Held–Suarez 多年气候态 | 多年积分达到统计平衡态未验证 |
| 14 | 长期积分守恒漂移监测 | 已做单步空间守恒验证（机器精度），未做长时间积分漂移监测 |

---

## 四、P8 Python 后处理收尾

| # | 工作项 | 位置 | 说明 |
|---|---|---|---|
| 15 | xarray+Dask 完整字段映射 | `python/cubedsphere/io.py` | 当前是骨架，未接真实 NetCDF 字段（`zonal_mean` 的维度名是假设的） |
| 16 | cartopy 地图投影图 | `python/cubedsphere/plotting.py` | 仅有 6 面板展开图，无地图投影（cartopy） |

---

## 已完成的核心成就（对照）

✅ **P0–P9 全部阶段框架**：模块化设计、接口先行、零警告编译、规范注释。

✅ **已完成的运行级/数值验证**（在当前环境完成的高价值部分）：
- 平流方案收敛阶实测：Center2=2.00 / Upwind3=2.98 / PPM=2.35 阶
- 平流通量配对守恒性实测：机器精度（1e-19 ~ 1e-14）
- 点积检验（TL/AD 转置）：rel = 2.65e-16
- 4D-Var 成本函数梯度（解析 vs 有限差分）：rel = 2.47e-10
- 全球面积积分收敛到 4π（相对误差随分辨率一阶收敛）
- CRC32C 校验（标准测试向量一致）
- 静力平衡、Kahan 求和、随机舍入、PPM 单调性等

---

## 建议优先级

**当前 Windows 环境可完成（无需 HPC）：**
1. 切线检验（#7）——用 `FiniteDifferenceTLM` 在纯本地验证手写 TL 的正确性
2. 背景误差递归滤波（#8）——纯本地可实现递归滤波算子
3. 后处理 I/O 闭环（#15、#16）——用 Zarr/纯 Python 侧先打通

**必须在 Linux HPC 完成（硬依赖）：**
- #1 完整构建、#2 MPI halo 交换、#3 NetCDF 写入、#4 Kokkos

**需完整时间积分跑通：**
- #11–#14 各标准算例与守恒漂移实测
