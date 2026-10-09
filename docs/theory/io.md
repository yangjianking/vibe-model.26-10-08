# 数据 I/O 子系统设计（P4 并行 I/O）

> 本文描述模式的输出、重启（检查点）与观测数据接口三大部分的设计依据与实现约定。
> 对应提示词集 P4 阶段：并行 NetCDF-4（HDF5 collective）、检查点/重启（N 进 M 出）、
> 观测数据接口（与 P5 同化联动）。

## 1. 总体职责划分

`src/io/` 下的 I/O 子系统遵循「接口先行」原则，分为三个正交模块：

| 模块 | 职责 | 后端 |
|---|---|---|
| `OutputWriter` / `NetCDFWriter` | 预报场输出 | NetCDF-4 并行；不可用时降级为二进制占位 |
| `RestartIO` | 检查点/重启 | 自描述二进制 + CRC32C 校验 |
| `ObsSpace` | 观测数据容器 + 质控 | 内存向量，供同化消费 |

三个模块互不依赖，便于独立替换后端（例如输出后端从 NetCDF 切换为 ADIOS2，不影响
重启与观测模块）。

## 2. 预报输出（OutputWriter）

### 2.1 输出约定

输出遵循 **CF Conventions 1.10**，网格元数据采用 **UGRID** 约定描述立方球拓扑：

- 变量命名优先使用 CF 标准名（`air_temperature`、`eastward_wind` 等）；
- 每个输出字段附带 `long_name` 与 `units`（见 `OutputField` 结构）；
- 网格元数据描述 6 个面板、面板间邻接关系、`ξ/η` 等角坐标与度量张量。

### 2.2 并行写出

NetCDF-4 底层为 HDF5，支持 `MPI_Collective` 的 `nc_put_vara_*` 并行写出。各 rank
仅持有本域（面板片段）的数据，通过 HDF5 hyperslab 选择写入全局网格对应区域：

```text
全局网格 = 6 面板 × (ncells × ncells × nlev)
   ↓ hyperslab 选择
rank k 写出本域 (panel, i0:i1, j0:j1, :) 片段
```

- **collective I/O**：所有 rank 协同写出，避免单点串行瓶颈；
- **压缩**：使用 HDF5 的 deflate / shuffle 过滤（参考 Klöwer et al. 2021 的
  位截断压缩思路，可在后续阶段加入有损压缩选项）。

### 2.3 降级路径

无 NetCDF 库时（如当前 Windows 开发环境），`NetCDFWriter` 退化为二进制占位输出
（`<path>.bin`），依次写 `(time, 字段尺寸, 字段数据)` 三元组。该路径仅用于调试与
接口联调，非生产格式。编译期通过 `CS_HAVE_NETCDF` 宏探测后端。

## 3. 重启 / 检查点（RestartIO）

### 3.1 设计目标

重启文件须满足**位级可复现**（bitwise reproducible），即从检查点恢复后继续积分，
结果与不中断的连续积分逐位一致。这要求：

1. **自描述元数据**：网格参数、坐标、halo、时间、迭代步、配置内容哈希全部落盘；
2. **固定序列化顺序**：8 个状态字段按固定顺序拼接为扁平缓冲；
3. **校验和**：CRC32C 覆盖状态字段，检测损坏。

### 3.2 文件格式（版本 1）

```text
┌─────────────────────────────────────────────┐
│ RestartMetadata 头部                        │
│   magic     = 0x43535253 ("CSRS")           │
│   version   = 1                             │
│   ncells / nhalo / nlev                     │
│   time / step                               │
│   config_hash（长度 + 内容）                │
├─────────────────────────────────────────────┤
│ 序列化状态字段（8 字段 × ncell 个 StateReal）│
│   ρ, ρu, ρv, ρw, ρθ, π', θ, Φ  (固定顺序)  │
├─────────────────────────────────────────────┤
│ CRC32C 校验和（覆盖状态字段）                │
└─────────────────────────────────────────────┘
```

### 3.3 CRC32C 校验

采用 **Castagnoli 多项式**（反射多项式 `0x82F63B78`），这是 CRC32C 的行业标准，
并可由 SSE4.2 的 `_mm_crc32_u64` 硬件指令加速（在 x86_64 HPC 节点上），软件逐位
实现作为通用回退。标准测试向量：`"123456789" → 0xE3069283`（已验证一致）。

### 3.4 N 进 M 出（rank 数可变重启）

重启文件以**全局网格**为存储粒度，与进程数无关。读取时各 rank 按当前域分解从全局
缓冲中切分自己所需的分片，因此支持 **N 进程写入 → M 进程读取**（N ≠ M）的弹性
重启。当前版本先实现单文件全局序列化，多文件并行分片在后续阶段补充。

## 4. 观测数据接口（ObsSpace）

### 4.1 观测抽象

`Observation` 结构统一描述一条观测：值、误差标准差 `σ_o`、位置（经/纬/高）、相对
同化窗起点的时间、观测类型与质控标记。首批观测类型见 `ObsType` 枚举：

| 枚举 | 含义 |
|---|---|
| `TEMP` | 探空温度 |
| `SYNOP` | 地面报 |
| `Satellite` | 卫星晴空辐射率 |
| `GNSS_RO` | GPS 无线电掩星折射率 |
| `Scatterometer` | 散射计海面风 |
| `AMV` | 大气运动矢量 |

### 4.2 时隙分箱

`time_slots(slot_width, window_len)` 将同化窗 `[0, window_len]` 划分为
`⌊window_len / slot_width⌋ + 1` 个时隙，把每条观测按其时间归入对应时隙。这是
4D-Var 在时间维度上组织观测向量 `y` 与观测算子 `H` 的基础。

### 4.3 质控（gross-error check）

`gross_error_check(k)` 做粗大误差检验：超出背景 `H(x_b) ± k·σ` 范围的观测标记为
拒绝（`qc_flag ≠ 0`）。当前版本用观测值自身的统计离群（均值 ± k·标准差）占位，
完整实现需接入背景场 `H(x_b)`（P5 同化阶段落地）。

## 5. 后续演进

- NetCDF-4 并行 collective 实际写入（HDF5 hyperslab + 压缩）；
- ADIOS2 后端（面向百 PB 级数据与异步 I/O）；
- 重启多文件并行分片与 N 进 M 出的完整实现；
- 有损压缩（Klöwer 2021 位截断）与位级可复现的权衡。

## 参考文献

- Eaton et al., *NetCDF Climate and Forecast (CF) Metadata Conventions*, v1.10.
- *UGRID Conventions* v1.0（非结构化网格元数据约定）。
- Klöwer, M. et al. (2021). *Compressing atmospheric data into its real information
  content.* Nat. Comput. Sci. **1**, 775–779.
- Castagnoli, G. et al. (1993). *Optimized CRC computation.* IEEE Trans. Comput.
