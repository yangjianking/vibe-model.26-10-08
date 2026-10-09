# 半隐式时间推进：离散化、稳定性与截断误差

> 模块：`src/timeint/` · 对应代码：`time_integrator.hpp/cpp`

## 1. 连续形式

将动力方程快波项与慢过程分离：

$$
\frac{\partial q}{\partial t} = \underbrace{\mathcal{L} q}_{\text{快波线性项}}
+ \underbrace{\mathcal{N}(q)}_{\text{慢非线性项}}
$$

## 2. 时间离散（Crank–Nicolson 隐式平均）

$$
\frac{q^{n+1} - q^n}{\Delta t} =
\frac{1}{2}\left[\mathcal{L} q^{n+1} + \mathcal{L} q^n\right]
+ \mathcal{N}(q^{n+1/2})
$$

即快波线性项用 Crank–Nicolson 平均（$\theta=1/2$），慢过程用显式或插值。
off-centering 参数 $\alpha\in[0.5,0.6]$ 通过 $\mathcal{L} q \to
\alpha\mathcal{L}q^{n+1} + (1-\alpha)\mathcal{L}q^n$ 引入，用于阻尼 $2\Delta t$
计算模态（见 §4）。

## 3. 截断误差阶数

- 快波线性项：Crank–Nicolson 平均为二阶精度 $O(\Delta t^2)$；
- 慢过程：取决于所配平流方案（中心二阶 / 三阶上游）。
- 半隐式消除声波 CFL 约束后，$\Delta t$ 仅受平流精度约束。

## 4. von Neumann 稳定性分析（概要）

对线性化方程做 von Neumann 分析，设 $\mathcal{L}$ 的特征值为纯虚数
$i\omega$（振荡模）。Crank–Nicolson 的放大因子：

$$
A = \frac{1 + i\alpha\omega\Delta t}{1 - i\alpha\omega\Delta t}
$$

对任意 $\omega$，$|A| = 1$，故半隐式对重力波/声波**无条件稳定**（中性稳定），
仅存 $2\Delta t$ 计算模态需 off-centering（$\alpha>0.5$）阻尼。平流 CFL 仅约束
半拉格朗日/有限体积平流的精度，而非稳定性。

## 5. 右端项构造

消元后 Helmholtz 方程的右端项由显式中间速度场的散度构造：

$$
R = -\beta\Delta t\, \nabla\cdot(\bar\rho\, \mathbf{v}^*)
$$

其中 $\mathbf{v}^*$ 为显式预报的中间速度（含平流、科氏力、垂直项，不含隐式
气压梯度），$\bar\rho$ 为参考态密度。散度用守恒形式离散（立方球度量，
`time_integrator.cpp::solve_helmholtz`），对应 `differential_operators.hpp` 的
通量散度。

## 6. 半隐式收益指标

- 显式 vs 半隐式最大稳定时间步之比 ≥ 4（提示词 P2 验收标准）；
- 10 天 RMS 差 < 1%（半隐式相对显式小步长参考解）。

## 7. 参考文献

- Robert (1981), *Atmos.-Ocean*, 19.
- Staniforth & Côté (1991), *Mon. Wea. Rev.*, 119.
- Wood et al. (2014), *Q. J. R. Meteorol. Soc.*, 140.
