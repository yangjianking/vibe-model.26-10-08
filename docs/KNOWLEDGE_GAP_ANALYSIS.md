# 补充知识与文档需求清单（含 ima 知识库可检索性评估）

> 目的：在项目附录 A「文献锚点清单」之外，系统梳理完成 P0–P9 各阶段**仍缺的  
> 知识与文档**，并逐项标注**能否通过 ima 知识库检索到**。  
> 生成日期：2026-10-08 · 审查者：AI-Agent（已实测检索 10 个知识库）

---

## 〇、结论速览

1. **ima 知识库覆盖面极广**，本项目附录 A 罗列的绝大多数文献（Sadourny、Ronchi、  
   Putman & Lin、Wood、Courtier、Giering、Váňa、Düben 等）均可在对应知识库检索到原文，  
   可作为「正文 → 原始出处」的核验闭环。
2. **缺口主要集中在三类**：
   - **工程 API 文档/手册**（Kokkos、NetCDF、ADIOS2、MPI、yaml-cpp、CMake、Spack）——  
     ima 知识库**基本不收录**，需走官方文档；
   - **特定数学公式/算法细节**（Helmholtz 谱分析、多重网格 LFA、切线性物理细节、  
     CVT 谱滤波、VarBC 系数推导）—— ima 可部分检索到论文，但**推导级细节需论文原文 +  
     专著**；
   - **可复现的软件工程规范**（位可重复归约、确定性并行）—— ima 覆盖较弱。

---

## 一、数学计算公式类（能否检索：部分可检索）

| 公式/算法                            | 项目位置                                 | ima 知识库可检索？                                           | 检索建议                              | 需额外提供                       |
| -------------------------------- | ------------------------------------ | ----------------------------------------------------- | --------------------------------- | --------------------------- |
| 立方球度量张量 $g_{ij},\sqrt G$ 完整推导    | `docs/theory/cubed_sphere_metric.md` | ✅ 可（Putman & Lin 2007、Ronchi 1996 原文）                 | 数值预报工具箱 KB 搜「cubed sphere metric」 | 已具备，无需                      |
| 协变/反变向量旋转与跨面板 halo               | `panel_rotation.hpp`                 | ✅ 可（Putman & Lin 2007 §3）                             | 同上                                | 已具备                         |
| Helmholtz 算子谱半径估计                | `docs/theory/helmholtz.md`           | ⚠️ 部分（Buckeridge 2010、GungHo）                         | GMD KB 搜「multigrid Helmholtz」     | 谱分析 + LFA 推导级细节             |
| 几何多重网格 LFA（局部 Fourier 分析）        | P2 收尾                                | ⚠️ 部分（Buckeridge & Scheichl 2010）                     | GMD/RMetS KB                      | LFA 教科书（Trottenberg 等 2001） |
| 半隐式 von Neumann 稳定性分析            | `docs/numerics/semi_implicit.md`     | ✅ 可（Staniforth & Côté 1991 综述）                        | RMetS KB                          | 已具备综述，缺显式放大因子推导             |
| PPM 单调限制器 + 截断误差阶                | `docs/numerics/advection.md`         | ✅ 可（Colella & Woodward 1984、LeVeque）                  | 数值预报工具箱 KB                        | 已具备                         |
| 背景误差协方差 B（NMC / 递归滤波 / 球谐谱）      | `src/da/`                            | ✅ 可（Parrish & Derber 1992、Bannister 2008、Fisher 2003） | ECMWF KB                          | Bannister 2008 综述原文         |
| CVT 控制变量变换（$v=B^{-1/2}\delta x$） | `four_dvar`                          | ⚠️ 部分                                                 | ECMWF KB 搜「control variable」      | Fisher 2003 完整推导            |
| 平衡算子（地转/统计回归）                    | `src/da/`                            | ⚠️ 部分                                                 | ECMWF KB                          | Derber & Bouttier 1999      |
| VarBC 变分偏差订正系数                   | `observation_operator`               | ✅ 可（Dee 2004、Auligné 2007）                            | ECMWF KB 搜「bias correction」       | Dee 2004 原文                 |
| 切线性物理（降水/边界层/扩散线性化）              | `src/da/tlm_ad`                      | ✅ 可（Janisková 1999/2013、Mahfouf 1999）                 | ECMWF KB                          | ECMWF TM 原文                 |
| 邻域法降水检验 FSS                      | P7                                   | ✅ 可（Roberts & Lean 2008）                              | RMetS KB                          | 已具备                         |
| 位涡拟能/角动量离散守恒（Thuburn 2008）       | `diagnostics`                        | ⚠️ 部分                                                 | GMD/RMetS KB                      | Thuburn 2008 原文             |

---

## 二、编程 API 文档 / 手册类（能否检索：**基本不可检索**，需官方文档）

ima 知识库以「学术文献 + 少量模式源码文档」为主，**几乎不含第三方库的官方 API  
手册**。以下全部需从官方站点获取：

| 库/工具                      | 项目用途               | 阶段    | ima 可检索？                         | 需提供的文档                                                                                      |
| ------------------------- | ------------------ | ----- | -------------------------------- | ------------------------------------------------------------------------------------------- |
| **Kokkos**                | 性能可移植内核抽象（P3 硬性要求） | P3    | ❌ 基本不可                           | Kokkos Core/Algorithm/Kernels 官方文档、Programming Guide、View/ExecPolicy/parallel_for/reduce 教程 |
| **NetCDF-C / NetCDF-C++** | 数据 I/O 主格式         | P4    | ❌ 不可                             | NetCDF-C API 手册、CF-1.10 约定文档、UGRID 元数据规范                                                    |
| **HDF5**                  | 底层并行 I/O           | P4    | ❌ 不可                             | HDF5 Parallel HDF5 教程（collective I/O）                                                       |
| **ADIOS2**                | 异步步进式输出            | P4    | ❌ 不可                             | ADIOS2 用户指南（BP5 引擎）                                                                         |
| **MPI (MPICH/OpenMPI)**   | 进程间通信、halo 交换      | P1/P4 | ❌ 不可                             | MPI-3.1 标准、MPI_Isend/Irecv 非阻塞通信最佳实践                                                        |
| **OpenMP**                | 节点内并行              | P3    | ❌ 不可                             | OpenMP 5.x 规范、SIMD/affinity/first-touch 章节                                                  |
| **yaml-cpp**              | 配置解析               | P9    | ❌ 不可                             | yaml-cpp 教程（官方 GitHub）                                                                      |
| **Catch2**                | 单元测试               | P9    | ❌ 不可                             | Catch2 官方文档                                                                                 |
| **CMake**                 | 构建系统               | P9    | ❌ 不可                             | CMake 官方手册（target/property/预设）、Profiling Guide                                              |
| **Spack**                 | HPC 工具链复现          | P9    | ❌ 不可                             | Spack 文档（spec 语法、environment）                                                               |
| **BLAS/LAPACK**           | 求解器底层              | P2    | ❌ 不可                             | LAPACK 用户指南                                                                                 |
| **Tapenade / Enzyme**     | TL/AD 自动微分交叉验证     | P6    | ❌ 不可                             | Tapenade 教程、Enzyme.jl 文档                                                                    |
| **RTTOV**                 | 卫星辐射传输             | P5    | ⚠️ 部分（工具箱有 rttov14_svr.pdf 用户手册） | RTTOV v14 用户手册                                                                              |
| **METIS/ParMETIS**        | 图划分                | P1    | ❌ 不可                             | METIS 手册                                                                                    |

> **结论**：Kokkos 编程（如你提到的）**确实需要单独提供 API 文档与手册**，ima  
> 知识库不收录此类第三方库官方文档。建议将 Kokkos Programming Guide、Kokkos  
> 官方 wiki、以及 Kokkos 在天气模式中的用法（LFRic/GungHo 的 PSyKAl）作为「工程  
> 参考文档」单独投喂，与学术文献（Edwards et al. 2014、Adams et al. 2019）区分开。



---

## 三、学术文献类（能否检索：**大部分可检索**）

对照提示词各阶段的文献锚点，ima 知识库的对应情况：

| 主题        | 附录 A 锚点                                                 | ima 对应知识库         | 可检索？  |
| --------- | ------------------------------------------------------- | ----------------- | ----- |
| 立方球/网格    | Sadourny 1972、Ronchi 1996、Putman & Lin 2007、Harris 2021 | 数值预报工具箱、RMetS     | ✅     |
| 半隐式/半拉格朗日 | Robert 1981、Tanguay 1990、Wood 2014                      | RMetS             | ✅     |
| 变分同化/伴随   | Courtier 1994、Giering 1998、Talagrand 1987               | ECMWF（海量）、RMetS   | ✅     |
| 混合精度      | Düben 2014、Váňa 2017、Nakano 2018                        | 数值预报工具箱、JGR/MWR   | ✅     |
| 检验算例      | Williamson 1992、Jablonowski 2006、Held-Suarez 1994、DCMIP | GMD、RMetS         | ✅     |
| Kokkos/众核 | Edwards 2014、Adams 2019                                 | GMD（gmd 相关篇）      | ✅（论文） |
| 教材专著      | Kalnay、Vallis、Coiffier、Jacobson 等                       | 数值预报工具箱（多部完整 PDF） | ✅     |

**需额外补充、ima 覆盖较弱或需精确到公式的**：

1. **Sadourny 1972 原文**（附录 A 已列但项目 references 未给卷期页码——已在本次增强中补全）。
2. **Thuburn 2008** 守恒离散综述（位涡拟能/角动量守恒）——GMD/RMetS 部分可检索。
3. **Janjić 1984** 地形斜率限制判据——ima 覆盖较弱，需原文。
4. **Klemp et al. 2008** 能量守恒形式 Rayleigh 海绵层——需原文。
5. **Davies 1976** 侧边界松弛——RMetS 可检索。
6. **Schmidt 1977** 网格拉伸——较老文献，ima 覆盖较弱。
7. **Arakawa 1966** 守恒离散经典、**Arakawa & Lamb 1977** C-grid——较老，需原文。

---

## 四、软件工程 / 规范类（能否检索：**基本不可**）

| 规范               | 用途          | ima 可检索？                              |
| ---------------- | ----------- | ------------------------------------- |
| CF-1.10 元数据约定    | NetCDF 属性规范 | ❌（需官方 CF 约定文档）                        |
| UGRID 网格元数据规范    | 立方球网格编码     | ❌（需官方 UGRID 文档）                       |
| BUFR 观测格式规范      | 观测读入        | ❌（需 WMO BUFR 手册；工具箱有 WMO-306 但需针对性章节） |
| 位可重复并行归约最佳实践     | 确定性并行       | ❌（需 Higham 2002 专著 + 实践文档）            |
| SLURM/PBS 作业脚本规范 | HPC 提交      | ❌（需集群管理员手册）                           |

---

## 五、总体建议：建议额外投喂的「知识包」清单

按优先级分三档：

### A 档（必须，影响能否继续编码）

1. **Kokkos 官方文档包**：Kokkos Core Programming Guide、Kokkos::View /  
   ExecutionPolicy / parallel_for / parallel_reduce / team 级 API、CUDA/HIP 后端配置。  
   —— ima 检索不到，需官方 wiki（kokkos.github.io）。
2. **并行 I/O 文档包**：NetCDF-C API、CF-1.10 + UGRID 约定、Parallel-HDF5 教程、  
   ADIOS2 BP5 指南。 —— ima 检索不到。
3. **MPI-3.1 + OpenMP 5.x 规范**（非阻塞通信、亲和性、SIMD）。 —— ima 检索不到。

### B 档（强烈建议，补公式推导细节）

1. **背景误差协方差 B 建模**：Bannister 2008（综述）、Fisher 2003（谱滤波）、  
   Derber & Bouttier 1999（平衡算子）。 —— ECMWF KB 可检索 Fisher/相关 TM。
2. **切线性物理线性化**：Janisková 1999/2013、Mahfouf 1999 的 ECMWF TM。 —— ECMWF KB 可检索。
3. **椭圆求解器与多重网格**：Buckeridge & Scheichl 2010 + Trottenberg 2001《Multigrid》。 —— 论文可检索，专著需补。
4. **守恒离散与诊断**：Thuburn 2008、Arakawa 1966。 —— 部分可检索，原文需补。
5. **切线性/伴随构建**：Giering & Kaminski 1998（完整 recipe）+ Griewank 2000（revolve）。 —— 可检索。

### C 档（建议，补工程与验证细节）

1. **VarBC 与质控**：Dee 2004、Andersson & Järvinen 1999、Desroziers 2005。 —— ECMWF KB 可检索。
2. **检验评估方法**：Jolliffe & Stephenson 2011（Forecast Verification 教材）、  
   Roberts & Lean 2008（FSS）。 —— 数值预报工具箱 KB 有完整教材 PDF。
3. **地形坐标与坡度限制**：Schär 2002、Janjić 1984、Klemp 2008。 —— 部分可检索。
4. **可复现与位一致**：Higham 2002 专著 + 混合精度实测（Váňa 2017、Zhang 2024）。 —— 可检索。

---

## 六、对「Kokkos 是否需要单独提供 API 文档」的直接回答

**需要，且是最高优先级。** 具体理由：

1. P3 明确「用 Kokkos 封装所有热点内核」为**硬性要求**，业务代码只允许用  
   `Kokkos::View / ExecPolicy / parallel_for / parallel_reduce`，这要求对 Kokkos  
   的执行空间（ExecutionSpace）、内存空间（MemorySpace）、内存布局（LayoutLeft/  
   LayoutRight）、team 级与 reduce 语义有精确掌握，而这些**不是学术论文能替代的**。
2. ima 知识库仅能检索到 Kokkos 的**学术文献**（Edwards et al. 2014 论文、LFRic 的  
   PSyKAl 用法论文），**检索不到 Kokkos 官方 API 手册与编程指南**。
3. 建议投喂的具体 Kokkos 文档：
   - Kokkos Core Programming Guide（View 模板、并行模式）
   - Kokkos 官方 wiki 的「Getting Started」与后端配置（CUDA/HIP/SYCL/OpenMP）
   - 天气/气候模式中的 Kokkos 用法案例：LFRic（Adams et al. 2019）、  
     GungHo 的 PSyKAl 分离思想（GMD KB 可检索论文）
