# 半隐式 Helmholtz 算子（Helmholtz Operator）

> 模块：`src/timeint/` · 对应代码：`elliptic_solver.hpp`、`bicgstab.cpp`

## 1. 问题陈述

半隐式时间推进将产生快波（重力波、声波）的线性项隐式处理，其余项显式。
对线性项做 Crank–Nicolson 时间平均后，可消元得到一个关于 Exner 压力扰动
$\delta\pi'$ 的三维椭圆方程（Helmholtz 方程），其高效求解是半隐式格式的
性能核心。

## 2. 从线性化方程到 Helmholtz 方程

### 2.1 线性化快波方程组

对非静力可压缩方程在参考态 $(\bar{\rho}, \bar{\theta}, \bar{\pi})$ 附近线性化，
快波耦合项为（气压梯度 ↔ 散度、浮力 ↔ 垂直速度）：

$$
\frac{\partial \delta\mathbf{v}}{\partial t} = -c_p\bar\theta\nabla(\delta\pi')
+ \dots
$$

$$
\frac{\partial \delta\pi'}{\partial t} = -\frac{c_s^2}{c_p\bar\theta}
\nabla\cdot(\bar\rho\,\delta\mathbf{v}) + \dots
$$

其中 $c_s$ 为参考态声速，$\bar\theta$ 为参考态位温。

### 2.2 Crank–Nicolson 时间平均与消元

对快波项用 Crank–Nicolson 隐式平均（权重 $\beta = \alpha\Delta t$，
$\alpha\in[0.5,0.6]$ 为 off-centering 参数），对 $\delta\pi'$ 做时间离散并
代入动量方程的隐式部分，消去 $\delta\mathbf{v}$，得到：

$$
\left( I - \beta^2 c_s^2\, \nabla\cdot\!\left(\bar\rho\, \nabla\right)\right)
\delta\pi' = R
$$

在等温参考态下简化为标准的 Helmholtz 方程：

$$
\left( I - \beta^2\Delta t^2 c_s^2\, \nabla^2 \right) \delta\pi' = R
$$

其中 $R$ 为显式残差项。

## 3. 立方球上的离散椭圆算子

在立方球 C 网格上，$\nabla^2$ 用度量张量离散（见 `cubed_sphere_metric.md` §5.3）：

$$
\nabla^2 \phi = \frac{1}{\sqrt{G}}\left[
\frac{\partial}{\partial\xi}\left(\sqrt{G}\,g^{11}\frac{\partial\phi}{\partial\xi}
+ \sqrt{G}\,g^{12}\frac{\partial\phi}{\partial\eta}\right)
+ \frac{\partial}{\partial\eta}\left(\sqrt{G}\,g^{21}\frac{\partial\phi}{\partial\xi}
+ \sqrt{G}\,g^{22}\frac{\partial\phi}{\partial\eta}\right)
\right]
$$

椭圆算子的各向异性系数（$g^{ij}$）随面板内位置变化，且跨面板时基向量旋转，
因此预条件器必须处理各向异性与面板耦合。

## 4. 求解策略

- **Krylov 子空间法**：BiCGStab（非对称算子）或 GMRES，matrix-free 实现，
  算子以 stencil 回调暴露（`EllipticSolver::ApplyOp`）。
- **预条件**：几何多重网格（各向异性系数处理）+ 垂直线松弛/ADI；
  周期方向可用 FFT。本阶段提供 Jacobi 预条件骨架与 BiCGStab 实现，
  多重网格作为可扩展策略预留。
- **可选**：JFNK（Jacobian-free Newton–Krylov，Knoll & Keyes 2004）用于
  强非线性场景。

## 5. 参考文献

- Robert, A. (1981). A stable numerical integration scheme for the primitive
  meteorological equations. *Atmos.-Ocean*, 19, 35–52.
- Tanguay, M., Robert, A., Laprise, R. (1990). A semi-implicit semi-Lagrangian
  fully compressible regional forecast model. *Mon. Wea. Rev.*, 118, 1970–1980.
- Wood, N., et al. (2014). *Q. J. R. Meteorol. Soc.*, 140, 1505–1528.
- Buckeridge, S., Scheichl, R. (2010). Parallel geometric multigrid for global
  weather prediction. *Q. J. R. Meteorol. Soc.*, 136, 2723–2735.
- Knoll, D. A., Keyes, D. E. (2004). Jacobian-free Newton-Krylov methods.
  *J. Comput. Phys.*, 193, 357–397.
