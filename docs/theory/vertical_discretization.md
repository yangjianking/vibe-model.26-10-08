# 垂直方向离散与静力平衡（Vertical Discretization）

> 模块：`src/dynamics/` · 对应代码：`equations.cpp::apply_vertical`

## 1. 问题陈述

非静力垂直动量方程需离散重力项与垂直气压梯度力。正确离散要求：在静力平衡
参考态下，垂直气压梯度力与重力**精确抵消**，避免产生虚假的垂直加速度。

## 2. 连续形式

垂直动量方程（非静力，守恒形式）：

$$
\frac{\partial(\rho w)}{\partial t} = -\frac{\partial p}{\partial z} - \rho g
+ \mathcal{N}_z
$$

其中 $\mathcal{N}_z$ 为垂直平流等非线性项。

## 3. Exner 压力形式的垂直气压梯度

用 Exner 压力 $\pi = (p/p_0)^\kappa$ 与位温 $\theta$ 表示：

$$
-\frac{\partial p}{\partial z} = -c_p\,\rho\,\theta\,\frac{\partial\pi}{\partial z}
- \rho g
$$

这一形式的动机（Harris et al. 2021 §4）：避免对密度的显式除法，且与热力学
方程的位温形式自洽，便于构造能量守恒配对。

## 4. 离散与静力平衡检验

垂直 Exner 压力梯度用相邻层中心差分（边界单侧差分），重力项 $-ρg$ 直接作用
于垂直动量。静力平衡要求：

$$
-c_p\,\bar\rho\,\bar\theta\,\frac{\partial\bar\pi}{\partial z} - \bar\rho g = 0
$$

已数值验证（`docs/theory/reference_state.md` 配套脚本）：等温参考态下合力
$\approx 10^{-15}$（机器精度），即参考态是垂直动量方程的精确稳态解。

## 5. 与半隐式的耦合

垂直气压梯度的隐式部分（$\partial\pi'/\partial z$ 与浮力项的耦合）由半隐式
Helmholtz 求解器处理（见 `helmholtz.md`）；显式部分在本模块计算。快波
（重力波、声波）的垂直传播通过隐式垂直算子稳定，避免 CFL 限制。

## 6. 参考文献

- Harris, L. M., et al. (2021). GFDL TM GFDL2021001, §5.
- Wood, N., et al. (2014). *Q. J. R. Meteorol. Soc.*, 140, 1505–1528.
