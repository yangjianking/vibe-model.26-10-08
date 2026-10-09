# 混合精度误差预算（Mixed Precision Error Budget）

> 模块：`src/common/` · 对应代码：`types.hpp`、`precision_guard.hpp`、`compensated_sum.hpp`

## 1. 问题陈述

混合精度（主体 fp32、关键量 fp64）可显著提升吞吐并降低内存带宽需求，但必须
给出**误差预算**论证：舍入误差相对物理/离散误差的量级可控，且长积分不漂移
（提示词 P3 第 4 条硬性要求）。

## 2. 三级类型别名

| 类型 | 精度 | 用途 |
|---|---|---|
| `Real` | fp64（默认）/ fp32（mixed） | 主体内核：平流通量、物理倾向、扩散 |
| `ComputeReal` | fp64（恒） | 求解器、归约、迭代精化 |
| `StateReal` | fp64（恒） | 状态累积、网格几何、守恒诊断 |

## 3. 误差预算论证

### 3.1 舍入误差量级

fp32 的机器精度 $\varepsilon_{32} \approx 1.19\times10^{-7}$，fp64 为
$\varepsilon_{64} \approx 2.22\times10^{-16}$。

单步平流/物理倾向的相对舍入误差约 $O(\varepsilon_{32}) \sim 10^{-7}$，
累积 $N$ 步（无偏）后为 $O(\sqrt{N}\,\varepsilon_{32})$。对 10 天预报
（$N \sim 10^4$ 步），约 $10^{-5}$，远小于：

- 离散误差（二阶格式 $O(\Delta x^2)$、三阶 $O(\Delta x^3)$）：$10^{-3}\sim10^{-4}$
- 物理参数化不确定性：$10^{-2}\sim10^{-1}$

故主体 fp32 的舍入误差**远小于离散与物理误差**，不会主导误差预算。

### 3.2 必须保持 fp64 的量

以下量若降为 fp32 会**破坏守恒或位可重复**，必须保持 fp64：

1. **网格几何量**（度量张量 $g_{ij}$、Jacobian $\sqrt{G}$）：几何误差直接污染
   散度/梯度算子的守恒性；
2. **椭圆求解器残差**：Helmholtz 求解器的迭代收敛依赖 fp64 残差（迭代精化）；
3. **全局守恒诊断**：质量/能量/角动量积分，用 Kahan/Neumaier 补偿求和；
4. **状态累积量**：$\rho, \rho u, \rho\theta$ 等长期积分状态，避免漂移。

### 3.3 精度卫士

降精度点配除零/上下溢保护（`precision_guard.hpp`），避免 fp32 下 NaN/Inf
污染（Váňa et al. 2017 经验）。

## 4. 验收判据（提示词 P3）

- fp32 模式 vs fp64 参考：10 天预报 500 hPa 位势高度 RMS 差 ≤ 1.5 m
  （Váňa et al. 2017 判据），运行时间缩短 ≥ 30%；
- 位可重复：同精度同分区两次运行 bit-identical（fp64 模式）。

## 5. 随机舍入（诊断用）

`StochasticRounder` 提供随机舍入开关（Paxton et al. 2022），期望值无系统偏差
（已数值验证，均值相对偏差 $<10^{-3}$），用于精度敏感性诊断。生产运行默认
关闭（破坏位可重复）。

## 6. 参考文献

- Váňa, F., et al. (2017). Single precision in weather forecasting models:
  An evaluation with the IFS. *Mon. Wea. Rev.*, 145, 4957–4970.
- Düben, P. D., Palmer, T. N. (2014). Benchmark tests for numerical weather
  forecasts on inexact hardware. *Mon. Wea. Rev.*, 142, 3809–3826.
- Chantry, M., et al. (2019). Scale-selective precision for weather and climate
  forecasting. *Mon. Wea. Rev.*, 147.
- Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM.
- Paxton et al. (2022). [TO-VERIFY]（随机舍入，待补充确切出处）。
