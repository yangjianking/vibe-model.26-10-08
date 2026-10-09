# 热力学关系与 Exner 压力（Thermodynamics）

> 模块：`src/dynamics/` · 对应代码：`thermo.hpp`

## 1. 干空气理想气体状态方程

$$
p = \rho R_d T
$$

其中 $R_d = 287.05\ \text{J}\cdot\text{kg}^{-1}\cdot\text{K}^{-1}$ 为干空气
气体常数。

## 2. 位温与 Exner 压力

位温定义：

$$
\theta = T\left(\frac{p_0}{p}\right)^{\kappa},
\qquad \kappa = \frac{R_d}{c_p} \approx 0.2857
$$

Exner 压力：

$$
\pi = \left(\frac{p}{p_0}\right)^{\kappa}
= \frac{T}{\theta}
$$

于是温度由位温与 Exner 压力恢复：$T = \theta\pi$。

## 3. 由预后变量恢复 $\pi$

预后变量为 $\rho$ 与 $\rho\theta$，则

$$
\theta = \frac{\rho\theta}{\rho}
$$

由状态方程与位温定义消去 $T$：

$$
p = \rho R_d T = \rho R_d \theta \pi
= \rho R_d \theta \left(\frac{p}{p_0}\right)^{\kappa}
$$

解得（Harris et al. 2021 §4）：

$$
\pi = \left(\frac{\rho\theta R_d}{\rho\, p_0}\right)
        ^{\kappa/(1-\kappa)}
$$

对应 `thermo.hpp::exner_from_rho_theta`。

## 4. 声速

$$
c_s^2 = \gamma R_d T,
\qquad \gamma = \frac{c_p}{c_v} \approx 1.4
$$

用于半隐式 Helmholtz 算子与 CFL 约束（见 `helmholtz.md`）。

## 5. 常数表（SI 单位）

| 符号 | 值 | 说明 |
|---|---|---|
| $R_d$ | 287.05 J/(kg·K) | 干空气气体常数 |
| $c_p$ | 1004.5 J/(kg·K) | 定压比热 |
| $c_v$ | 717.5 J/(kg·K) | 定容比热 |
| $\kappa$ | 0.2857 | $R_d/c_p$ |
| $p_0$ | 100000 Pa | 参考气压 |
| $g$ | 9.80665 m/s² | 重力加速度 |

## 6. 参考文献

- Harris, L. M., et al. (2021). GFDL TM GFDL2021001, §4.
- Dutton, J. A. (1986). *The Ceaseless Wind: An Introduction to the Theory of
  Atmospheric Motion*. Dover.
