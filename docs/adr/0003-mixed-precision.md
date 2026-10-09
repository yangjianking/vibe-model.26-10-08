# ADR-0003：混合精度策略 —— 三级类型别名 + 补偿求和

- **状态**：已采纳
- **日期**：2026-10-08

## 背景

提示词 P3 要求混合精度抽象：主体计算可 fp32，网格几何/求解器/守恒诊断保持
fp64，且必须给出误差预算论证。

## 决策

采用三级类型别名 `Real / ComputeReal / StateReal`：
- `Real`：主体内核（平流/物理/扩散），CMake 选项 `CS_PRECISION` 切换 fp64/mixed；
- `ComputeReal` 与 `StateReal`：恒 fp64，用于求解器、网格几何、状态累积、
  守恒诊断。
- 全局归约用 Kahan/Neumaier 补偿求和（禁朴素求和）。
- 降精度点配精度卫士（除零/上下溢保护），随机舍入作诊断开关（默认关闭）。

## 后果

- 误差预算论证见 `docs/numerics/mixed_precision.md`：fp32 舍入误差
  $O(\sqrt{N}\varepsilon_{32})\sim10^{-5}$，远小于离散误差 $10^{-3}$ 与物理
  不确定性 $10^{-2}$。
- 位可重复：fp64 模式 + 确定性归约树保证；mixed 模式不保证位可重复（预期）。

## 参考文献

- Váňa et al. (2017), Düben & Palmer (2014), Chantry et al. (2019), Higham (2002)。
