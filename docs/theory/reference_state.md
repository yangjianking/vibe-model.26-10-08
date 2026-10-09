# 参考态与半隐式线性化（Reference State）

> 模块：`src/dynamics/` · 对应代码：`reference_state.hpp/cpp`

## 1. 问题陈述

半隐式时间推进需将非静力可压方程在参考态附近线性化，把产生快波的线性项
（气压梯度 ↔ 散度耦合、浮力 ↔ 垂直速度耦合）隐式处理。参考态取静力平衡
大气，预后变量分解为参考态 + 扰动：

$$
\rho = \bar\rho + \rho',\qquad
\rho\theta = \bar\rho\bar\theta + (\rho\theta)',\qquad
\pi = \bar\pi + \pi'
$$

## 2. 等温静力平衡参考态

取等温参考温度 $T = \bar T$（常数），静力平衡

$$
\frac{\partial \bar p}{\partial z} = -\bar\rho g,\qquad
\bar p = \bar\rho R_d \bar T
$$

解析解为指数剖面：

$$
\bar p(z) = p_s e^{-z/H},\qquad
\bar\rho(z) = \rho_s e^{-z/H},\qquad
H = \frac{R_d \bar T}{g}
$$

其中 $H$ 为标高。对 $\bar T = 250$ K，$H \approx 7317$ m。

## 3. 参考态派生量

- Exner 压力：$\bar\pi = (\bar p/p_0)^{\kappa}$
- 位温：$\bar\theta = \bar T / \bar\pi$（等温大气位温随高度递增，静力稳定）
- 声速平方：$c_s^2 = \gamma R_d \bar T$（等温下为常数，$\gamma = c_p/c_v$）

## 4. 半隐式线性化与 Helmholtz 右端项

将快波线性项用 Crank–Nicolson 隐式平均后，对 $\pi'$ 消元得到：

$$
\left(I - \beta^2\Delta t^2\, c_s^2\, \nabla^2\right) \pi' = R
$$

右端项 $R$ 由显式残差中的散度项构造：

$$
R = -\beta\Delta t\,\nabla\cdot(\bar\rho\, \mathbf{v}^*) + \dots
$$

其中 $\mathbf{v}^*$ 为显式预报的中间速度。完整对接见
`docs/numerics/semi_implicit.md`（当前 `time_integrator.cpp::solve_helmholtz`
已用参考态声速构造算子，右端项 $R$ 与动力散度的精确对接为后续阶段收尾项）。

## 5. 参考文献

- Harris, L. M., et al. (2021). GFDL TM GFDL2021001, §4–5.
- Wood, N., et al. (2014). *Q. J. R. Meteorol. Soc.*, 140, 1505–1528.
- Staniforth, A., Côté, J. (1991). *Mon. Wea. Rev.*, 119, 2206–2223.
