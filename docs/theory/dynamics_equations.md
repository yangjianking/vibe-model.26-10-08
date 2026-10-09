# 非静力动力方程与守恒配对（Dynamics Equations）

> 模块：`src/dynamics/` · 对应代码：`equations.hpp`、`equations.cpp`

## 1. 问题陈述

实现非静力可压缩大气动力核心。静力近似作为编译/配置开关保留，本框架默认
非静力。预后变量采用密度加权形式，以保证质量、动量、能量的通量守恒离散。

## 2. 连续方程与量纲

预后变量（Harris et al. 2021 §3）：

| 变量 | 含义 | 量纲 |
|---|---|---|
| $\rho$ | 干空气质量密度 | kg·m⁻³ |
| $\rho u$ | ξ 方向动量密度（协变） | kg·m⁻²·s⁻¹ |
| $\rho v$ | η 方向动量密度（协变） | kg·m⁻²·s⁻¹ |
| $\rho w$ | 垂直动量密度 | kg·m⁻²·s⁻¹ |
| $\rho\theta$ | 位温密度 | K·kg·m⁻³ |

### 2.1 质量守恒（连续方程，通量形式）

$$
\frac{\partial \rho}{\partial t} = -\nabla\cdot(\rho \mathbf{v})
$$

### 2.2 动量方程（守恒形式）

$$
\frac{\partial (\rho \mathbf{v})}{\partial t} =
-\nabla\cdot(\rho \mathbf{v}\otimes\mathbf{v})
- \nabla p - \rho g\,\hat{k}
- 2\rho\,\boldsymbol{\Omega}\times\mathbf{v}
+ \mathcal{F}_{\text{phys}}
$$

其中 $\boldsymbol{\Omega}$ 为地球自转角速度向量，$\mathcal{F}_{\text{phys}}$
为物理过程（摩擦、重力波拖曳等）源项。

### 2.3 热力学方程（位温密度守恒）

$$
\frac{\partial (\rho\theta)}{\partial t} =
-\nabla\cdot(\rho\theta \mathbf{v}) + \mathcal{Q}_{\text{phys}}
$$

$\mathcal{Q}_{\text{phys}}$ 为非绝热加热（辐射、潜热释放等）。

## 3. 坐标与离散化

水平：立方球 gnomonic equiangular 投影（见 `cubed_sphere_metric.md`），
Arakawa C-grid 交错。垂直：地形跟随坐标（可选纯高度坐标）。

C 网格交错约定：

| 量 | 位置 |
|---|---|
| 标量 $\rho,\theta,\pi$ | 单元中心 $(i,j)$ |
| $u$（ξ 速度） | 东面 $(i+1/2, j)$ |
| $v$（η 速度） | 北面 $(i, j+1/2)$ |
| 涡度/散度 | 顶点 $(i+1/2, j+1/2)$ |

## 4. 通量与守恒配对

### 4.1 气压梯度与散度的能量配对（Arakawa–Lamb 1981）

离散总能量守恒要求：气压梯度力（作用于面）与散度算子（作用于单元中心）
采用**互为伴随**的离散形式。在立方球 C 网格上，这等价于：

- 散度用守恒通量形式（§5.1 度量文档）；
- 气压梯度力用与散度严格配对的面差分形式（见
  `docs/numerics/pressure_gradient.md`）。

气压梯度采用 Exner 压力 $\pi$ 与位温 $\theta$ 的耦合形式：

$$
-\frac{1}{\rho}\nabla p \;\longrightarrow\; -c_p\theta\,\nabla\pi
$$

其中 $\pi = (p/p_0)^{\kappa}$。这避免了对密度的显式除法，且与热力学方程
的位温形式自洽（Harris et al. 2021 §4）。

### 4.2 离散总能量守恒的代数证明（概要）

对动量方程点乘 $\rho\mathbf{v}$ 并对全球求和，对热力学方程乘 $c_p\theta$
求和，两式相加，利用散度/梯度的伴随配对，边界通量项经面板连续性与
halo 交换抵消，得到：

$$
\frac{d}{dt}\int \left(\frac{1}{2}\rho|\mathbf{v}|^2 + c_v\rho T +
\rho gz\right)dV = 0
$$

完整逐项证明见 `docs/numerics/pressure_gradient.md`（后续阶段补全代数细节）。

## 5. 科氏力与曲率项

球面上协变形式的科氏力包含度量 Christoffel 项（Putman & Lin 2007 §3）。
本框架在 `apply_coriolis` 中给出 $\beta$-平面近似骨架（$f = 2\Omega\sin\phi$），
完整曲率项在后续阶段补充。

## 6. 参考文献

- Arakawa, A., Lamb, V. R. (1981). A potential enstrophy and energy conserving
  scheme for the shallow water equations. *Mon. Wea. Rev.*, 109, 18–36.
- Lin, S.-J., Rood, R. B. (1996). Multidimensional flux-form semi-Lagrangian
  transport schemes. *Mon. Wea. Rev.*, 124, 2046–2070.
- Harris, L. M., et al. (2021). GFDL TM GFDL2021001.
- Wood, N., et al. (2014). An inherently mass-conserving semi-implicit
  semi-Lagrangian discretization of the deep-atmosphere global non-hydrostatic
  equations. *Q. J. R. Meteorol. Soc.*, 140, 1505–1528.
