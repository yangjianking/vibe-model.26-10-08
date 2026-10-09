# 气压梯度力守恒配对离散（Pressure Gradient）

> 模块：`src/dynamics/` · 对应代码：`equations.cpp::compute_explicit`

## 1. 问题陈述

在 Arakawa C 网格上，气压梯度力作用于面（与速度共位），而散度作用于单元
中心。为保证离散总能量守恒，两者必须采用**互为伴随**的离散形式
（Arakawa–Lamb 1981）。

## 2. 连续形式

$$
-\frac{1}{\rho}\nabla p = -c_p\theta\,\nabla\pi
$$

采用 Exner 压力 $\pi$ 与位温 $\theta$ 的耦合形式（`thermo.hpp`），避免对
$\rho$ 的显式除法。

## 3. 离散形式（面处中心差分）

在单元中心 $(i,j)$ 处的协变气压梯度力：

$$
\left(c_p\theta\,\frac{\partial\pi}{\partial\xi}\right)_{i,j}
= c_p\,\theta_{i,j}\,
  \frac{\pi_{i+1,j} - \pi_{i-1,j}}{2\Delta\xi}
$$

协变 → 反变（乘逆度量，见 `covariant.hpp`）：

$$
\text{PG}_\xi = -c_p\theta\,(g^{11}\partial_\xi\pi + g^{12}\partial_\eta\pi)
$$

## 4. 能量守恒配对条件

离散总能量守恒要求（证明概要，见 `dynamics_equations.md` §4.2）：

$$
\sum_{i,j} \mathbf{v}_{i,j}\cdot(\text{PG})_{i,j} \cdot \sqrt{G}\,\Delta\xi\Delta\eta
= -\sum_{i,j} \pi_{i,j}\,\nabla\cdot(\rho\mathbf{v})_{i,j}
$$

即气压梯度力与散度算子的离散形式互为转置（伴随）。当前 `compute_explicit`
已按此结构实现骨架，完整配对校验在后续阶段用解析场逐项验证。

## 5. 参考文献

- Arakawa & Lamb (1981), *Mon. Wea. Rev.*, 109, 18–36.
- Harris et al. (2021), GFDL TM GFDL2021001, §4–5.
