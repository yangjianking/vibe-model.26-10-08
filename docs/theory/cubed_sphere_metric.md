# 立方球网格度量张量推导（Cubed-Sphere Metric）

> 模块：`src/grid/` · 对应代码：`include/cubed_sph/grid/metric.hpp`、`covariant.hpp`
> 本文档给出 gnomonic equiangular 立方球投影的完整度量推导，是散度/梯度/涡度
> 差分算子（`differential_operators.hpp`）的几何基础。

## 1. 问题陈述

全球大气模式需要在球面上离散水平方程。经典经纬网格存在极点奇异（极区经线汇聚），
导致时间步长受 CFL 严格限制。立方球网格将球面划分为 6 个近似等面积的方形面板，
每个面板用局部笛卡尔坐标描述，从根本上消除极点奇异，并天然支持全局无缝网格
（Sadourny 1972；Ronchi et al. 1996）。

## 2. gnomonic equiangular 投影

### 2.1 定义

设面板中心方向为单位向量 $\hat{n}$（面板法向），面板局部正交基为
$(\hat{e}_1, \hat{e}_2, \hat{n})$。面板局部**等角坐标** $(\xi, \eta)$ 定义为

$$
\xi = \tan\theta_\xi,\qquad \eta = \tan\theta_\eta,
\qquad \theta_\xi, \theta_\eta \in [-\pi/4, \pi/4]
$$

因此 $\xi, \eta \in [-1, 1]$（**无量纲**，面板角点对应 $\xi = \eta = \pm 1$）。
球面点的映射（gnomonic 投影）为：

$$
\mathbf{x}(\xi, \eta) = \frac{\hat{n} + \xi\,\hat{e}_1 + \eta\,\hat{e}_2}{r},
\qquad r = \sqrt{1 + \xi^2 + \eta^2}
$$

其中 $\mathbf{x}$ 是球面单位向量（半径已归一化到 1）。当 $\xi = \eta = 0$ 时
$\mathbf{x} = \hat{n}$，即面板中心。该投影把立方体面中心投影到球面，四个角点
对应 $\xi, \eta = \pm 1$（即 $\theta = \pm \pi/4$）。

### 2.2 等面积性质

投影的 Jacobian（面积元缩放）为

$$
\sqrt{G} = \left| \frac{\partial \mathbf{x}}{\partial \xi} \times
\frac{\partial \mathbf{x}}{\partial \eta} \right| = \frac{1}{r^3}
$$

其中 $r^2 = 1 + \xi^2 + \eta^2$。因此面积元

$$
dA = \sqrt{G}\; d\xi\, d\eta = \frac{d\xi\, d\eta}{(1+\xi^2+\eta^2)^{3/2}}
$$

在面板中心（$r=1$）面积元最大，向面板边缘（$r\to\sqrt{3}\approx 1.732$）
递减，6 个面板合计精确覆盖整个球面（$4\pi$ 球面度）。相比等距投影，equiangular
投影的面板内网格更均匀，面积变形更小。

> **验证**：$\iint_{\text{面板}} \sqrt{G}\,d\xi d\eta =
> \int_{-1}^{1}\int_{-1}^{1}(1+\xi^2+\eta^2)^{-3/2}d\xi d\eta = 2\pi/3$，
> 6 面板合计 $4\pi$（数值积分相对误差 $< 10^{-6}$，见
> `tests/unit/test_grid_metric.cpp` 的配套验证脚本）。

## 3. 度量张量

### 3.1 协变度量分量 $g_{ij}$

度量张量由切向量内积定义：

$$
g_{11} = \frac{\partial \mathbf{x}}{\partial \xi} \cdot
         \frac{\partial \mathbf{x}}{\partial \xi},
\qquad
g_{22} = \frac{\partial \mathbf{x}}{\partial \eta} \cdot
         \frac{\partial \mathbf{x}}{\partial \eta},
\qquad
g_{12} = \frac{\partial \mathbf{x}}{\partial \xi} \cdot
         \frac{\partial \mathbf{x}}{\partial \eta}
$$

对投影式求导（利用 $r = \sqrt{1+\xi^2+\eta^2}$，$\partial r/\partial\xi = \xi/r$）：

$$
\frac{\partial \mathbf{x}}{\partial \xi}
= \frac{\hat{e}_1}{r} - \frac{(\hat{n}+\xi\hat{e}_1+\eta\hat{e}_2)\xi}{r^3}
$$

代入并利用基向量正交性（$\hat{n}\cdot\hat{e}_1=0$ 等）与归一化（$\hat{e}_i\cdot\hat{e}_i=1$），
得到（Putman & Lin 2007, Eq. A1–A6；Ronchi et al. 1996, Eq. 7–10）：

$$
g_{11} = \frac{1 + \eta^2}{r^4},
\qquad
g_{22} = \frac{1 + \xi^2}{r^4},
\qquad
g_{12} = -\frac{\xi\eta}{r^4}
$$

### 3.2 行列式与 Jacobian

$$
\det(g) = g_{11}g_{22} - g_{12}^2
= \frac{(1+\eta^2)(1+\xi^2) - \xi^2\eta^2}{r^8}
= \frac{1+\xi^2+\eta^2}{r^8}
= \frac{1}{r^6}
$$

因此 $\sqrt{G} = \sqrt{\det(g)} = 1/r^3$，与 §2.2 的几何推导一致。

### 3.3 逆度量张量 $g^{ij}$

$$
g^{11} = \frac{g_{22}}{\det(g)} = (1+\xi^2)\, r^2,
\qquad
g^{22} = \frac{g_{11}}{\det(g)} = (1+\eta^2)\, r^2,
\qquad
g^{12} = -\frac{g_{12}}{\det(g)} = \xi\eta\, r^2
$$

## 4. 协变 / 反变分量变换

切平面内任意向量 $\mathbf{u}$ 可表示为

$$
\mathbf{u} = u^{\xi}\,\partial_\xi + u^{\eta}\,\partial_\eta
            = u_\xi\, g^{\xi j}\partial_j + \dots
$$

其中 $u^i$ 为**反变分量**，$u_i$ 为**协变分量**，满足

$$
u^i = g^{ij} u_j,
\qquad
u_i = g_{ij} u^j
$$

对应代码 `covariant.hpp` 中的 `covariant_to_contravariant` /
`contravariant_to_covariant`。

**采用协变分量存储的动机**（Putman & Lin 2007 §3）：

1. 散度算子 $\nabla\cdot\mathbf{u}$ 在协变分量下形式简洁且严格守恒；
2. 跨面板时协变分量按局部基向量旋转，旋转矩阵仅依赖面板朝向，与位置无关，
   极大简化 halo 交换（见 `panel_rotation.hpp`）。

## 5. 差分算子中的度量应用

### 5.1 散度（守恒形式）

$$
\nabla \cdot \mathbf{u} = \frac{1}{\sqrt{G}}\left[
  \frac{\partial}{\partial \xi}\left(\sqrt{G}\, u^\xi\right) +
  \frac{\partial}{\partial \eta}\left(\sqrt{G}\, u^\eta\right)
\right]
$$

对应 `differential_operators.hpp::divergence`，其中 $\sqrt{G}\,u^\xi$ 为
"通量密度"（`covariant.hpp::flux_density`），是协变→反变后乘 $\sqrt{G}$ 得到。

### 5.2 涡度（协变形式）

$$
\zeta = \frac{1}{\sqrt{G}}\left(
  \frac{\partial v_\xi}{\partial \xi} - \frac{\partial u_\eta}{\partial \eta}
\right)
$$

用协变分量直接计算，避免了 Christoffel 符号的显式处理
（Putman & Lin 2007）。

### 5.3 拉普拉斯算子

$$
\nabla^2 \phi = \frac{1}{\sqrt{G}}\left[
  \frac{\partial}{\partial \xi}\left(\sqrt{G}\, g^{11}
  \frac{\partial \phi}{\partial \xi} + \sqrt{G}\, g^{12}
  \frac{\partial \phi}{\partial \eta}\right) + \dots
\right]
$$

用于半隐式 Helmholtz 算子与散度阻尼（见 `docs/theory/helmholtz.md`）。

## 6. 面板朝向与向量旋转

6 个面板的面板中心方向约定（与 Harris et al. 2021 一致）：

| 面板 | 中心方向 | $\hat{e}_1$ | $\hat{e}_2$ |
|---|---|---|---|
| Face1 | $+x$ | $+y$ | $+z$ |
| Face2 | $-y$ | $-x$ | $+z$ |
| Face3 | $-x$ | $-y$ | $+z$ |
| Face4 | $+y$ | $+x$ | $+z$ |
| North | $+z$ | $+x$ | $+y$ |
| South | $-z$ | $+x$ | $-y$ |

跨面板传递向量时，旋转矩阵 $R$ 的元素为两面板基向量的内积：

$$
R = \begin{pmatrix}
\hat{e}_1^{(s)}\cdot\hat{e}_1^{(d)} & \hat{e}_2^{(s)}\cdot\hat{e}_1^{(d)} \\
\hat{e}_1^{(s)}\cdot\hat{e}_2^{(d)} & \hat{e}_2^{(s)}\cdot\hat{e}_2^{(d)}
\end{pmatrix}
$$

对应 `panel_rotation.hpp::panel_rotation`。

## 7. 实现要点与代码映射

| 概念 | 代码位置 |
|---|---|
| gnomonic 投影 | `src/grid/cubed_sphere_grid.cpp::compute_metric` |
| 度量缓存 | `CubedSphereGrid::build_metrics` |
| 协变/反变变换 | `grid/covariant.hpp` |
| 面板旋转 | `src/grid/panel_rotation.cpp` |
| 差分算子 | `grid/differential_operators.hpp` |
| 面板拓扑 | `src/grid/cubed_sphere_grid.cpp::init_neighbors` |

## 8. 参考文献

- Sadourny, R. (1972). Conservative finite-difference approximations of the
  primitive equations on quasi-uniform spherical grids. *Mon. Wea. Rev.*, 100, 136–144.
- Ronchi, C., Iacono, R., Paolucci, P. S. (1996). The "cubed sphere": a new method
  for the solution of partial differential equations in spherical geometry.
  *J. Comput. Phys.*, 124, 93–114.
- Putman, W. M., Lin, S.-J. (2007). Finite-volume transport on various cubed-sphere
  grids. *J. Comput. Phys.*, 227, 55–78.
- Harris, L. M., et al. (2021). A Scientific Description of the GFDL FV³ Dynamical
  Core. *GFDL Tech. Memo.* GFDL2021001, §3.
