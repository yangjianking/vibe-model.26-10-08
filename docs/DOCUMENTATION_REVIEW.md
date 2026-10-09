# 文档内容增强分析报告（基于 ima 知识库 + 提示词集的对照审查）

> 审查对象：`vibe-model.2026-10-08` 立方球 NWP 模式项目
> 审查依据：《NWP 模式开发 · AI-Agent 提示词集 v2.0》P0–P9 各阶段验收标准、
> 附录 A 文献锚点清单，以及 ima 知识库可检索到的权威文献。
> 生成日期：2026-10-08

---

## 一、审查结论总览

现有文档体系（21 个文件，`docs/theory` + `docs/numerics` + `docs/adr` +
`docs/user`）**框架完整、科学严谨**，已覆盖 P1–P2 核心推导，并具备 P3–P6 的
接口层文档。对照提示词 P9 第 4 条「每个算法文档固定小节」的规范，现有文档在
**数学推导链完整性、文献可查证性、跨知识库关联**三方面仍有可增强空间。

| 维度 | 现状评估 | 增强优先级 |
|---|---|---|
| 公式推导（连续→离散→截断误差→稳定性） | 部分文档缺「截断误差阶数」与「von Neumann 分析」完整小节 | 高 |
| 文献可查证性 | 少数 [TO-VERIFY] 未闭环，缺 DOI | 中 |
| 与 ima 知识库文献映射 | 未显式标注「可检索文献」 | 中 |
| 验收标准与测试对应 | 分散在各文档，缺统一追溯表 | 中 |

---

## 二、逐文档增强建议（含知识库可检索到的对应文献）

### 1. `docs/numerics/advection.md` — 补充截断误差阶数推导与文献

当前已给出 PPM 重构，但可补：
- **中心二阶 / 三阶上游 / PPM 的截断误差阶数**显式推导（modifed equation），
  与 `python/cubedsphere/convergence.py` 实测（Center2=2.0 / Upwind3=3.0 /
  PPM=2.35）一一对应。
- **正性与单调性**：PPM 限制器（Colella & Woodward 1984）的 TVD/正性证明。

**ima 可检索文献**（数值预报工具箱 KB `7295373076860673`）：
- `Finite Volume Methods for Hyperbolic Problems - LeVeque`（有限体积经典教材）
- `Numerical Methods for Wave Equations in Geophysical Fluid - Durran 1999`
- `Semi-Lagrangian Advection Methods... - Steven James Fletcher 2019`
- Huang-2022-QJ（多矩有限体积动力核）、Chen-2008/2014/2023-JCP、Tang 系列

### 2. `docs/numerics/semi_implicit.md` — 补充 von Neumann 稳定性分析

提示词 P2 第 3 条要求「给出 von Neumann 稳定性分析（半隐式对重力波无条件稳定）」。
建议补：
- 线性化快波系统（重力波 + 声波）的半隐式 Crank–Nicolson 放大因子推导；
- off-centering 参数 $\alpha \in [0.5, 0.6]$ 对 $2\Delta t$ 波的数值阻尼显式表达。

**ima 可检索文献**（RMetS KB `7345834534381192` + GMD KB `7311002446073452`）：
- Wood et al. 2014 (ENDGame, QJRMS 140)
- Staniforth & Côté 1991 综述、Côté et al. 1998 (GEM)

### 3. `docs/theory/reference_state.md` / `helmholtz.md` — 谱性质与预条件

- 补充离散 Helmholtz 算子 $I - \beta^2\Delta t^2 c_s^2 \nabla^2$ 的谱半径估计
  与 BiCGStab 收敛因子；多重网格预条件的 LFA（局部 Fourier 分析）框架。

**ima 可检索文献**（GMD KB）：
- Buckeridge & Scheichl 2010（并行几何多重网格，QJRMS）
- GungHo 系列（gmd-15-6601 等）、Dedner et al. 2016

### 4. `docs/theory/four_dvar.md` / `tlm_ad.md` — 补充闭环文献

当前点积检验、梯度检验已通过，可补：
- 增量 4D-Var 内/外循环与多增量（multi-incremental）的收敛性论证；
- 伴随「先离散后线性化」原则与 checkpointing（revolve）的完整引文。

**ima 可检索文献**（ECMWF KB `7296207936294738` 极为丰富）：
- Courtier et al. 1994、Rabier et al. 2000、Veersé & Thépaut 1998
- Giering & Kaminski 1998（伴随代码构建）、Griewank & Walther 2000（revolve）
- 大量 ECMWF 4D-Var/TL/AD 技术备忘录（incremental、long-window、convergence 等）

### 5. `docs/numerics/mixed_precision.md` — 补误差预算实测锚点

- 补 Zhang et al. 2024 (GRIST 混合精度, GMD 17) 作为国产模式混合精度实测锚点；
- 补 Düben & Palmer 2014、Nakano 2018、Chantry 2019 的完整引用。

**ima 可检索文献**（数值预报工具箱 + GMD KB）均有收录。

### 6. 新增 `docs/theory/physics_parameterizations.md`（建议）

当前物理参数化为 Dummy 占位（TODO #5），提示词 P2 要求物理接口可扩展。
建议新增物理参数化设计文档，覆盖辐射/边界层/微物理/重力波拖曳接口。

**ima 可检索文献**（数值预报工具箱 KB 含 GFS/FV3 物理源码文档）：
- `gfs_phys_time_vary.md`、`PBL.md`、`MYNN_EDMF.txt`、`UGWPv0.txt`、
  `GFDL_cloud.txt`、`rttov14_svr.pdf`（辐射传输）等。

---

## 三、建议新增文档清单

| # | 建议文档 | 对应阶段 | 必要性 |
|---|---|---|---|
| 1 | `docs/theory/physics_parameterizations.md` | P2/P3 | 物理接口设计依据 |
| 2 | `docs/numerics/stability_analysis.md` | P2 | 集中放 von Neumann 分析 |
| 3 | `docs/theory/background_error_covariance.md` | P5 | B 矩阵（NMC/集合）设计 |
| 4 | `docs/adr/0004-parallel-io-backend.md` | P4 | ADIOS2 vs NetCDF-4 决策 |
| 5 | `docs/numerics/multigrid.md` | P2 收尾 | 多重网格预条件 LFA |
| 6 | `docs/theory/halo_exchange.md` | P1 收尾 | 跨面板向量旋转与守恒 |

---

## 四、验收标准 ↔ 文档追溯表（建议维护）

提示词 P9 要求「文档中每个离散公式编号与代码注解一一对应」，建议在 README
或独立 `docs/traceability.md` 建立「验收标准 → 文档 → 代码 → 测试」四列追溯表，
一键核查每个验收点的闭环状态（当前分散在 TODO.md 与各文档）。
