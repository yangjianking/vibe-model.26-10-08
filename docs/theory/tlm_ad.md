# 切线性模型与伴随模型（P6 手写 TL/AD）

> 本文描述切线性模型（TL）与伴随模型（AD）的数学定义、手写实现策略与
> 点积检验。对应提示词集 P6 阶段。

## 1. 切线性模型（Tangent Linear, TL）

预报模型 $M: x(t_0) \mapsto x(t_1)$ 的切线性模型是它对初始条件的 Gateaux
导数（线性化）：

$$
M'(x) \, \delta x = \lim_{\varepsilon \to 0}
\frac{M(x + \varepsilon\,\delta x) - M(x)}{\varepsilon}
$$

对显式时间推进 $M(x) = x + \Delta t\, F(x)$，切线性传播为：

$$
\delta x \leftarrow \delta x + \Delta t\, F'(x)\, \delta x
$$

其中 $F'(x)$ 是动力倾向 $F$ 的 Jacobian。

## 2. 伴随模型（Adjoint, AD）

伴随模型是切线性算子的转置（对欧氏内积 $\langle\cdot,\cdot\rangle$）：

$$
\langle M'(x)\,u,\; v\rangle = \langle u,\; M'(x)^\mathsf{T}\, v\rangle
$$

伴随传播沿时间反向，用于 4D-Var 的梯度反推。对显式推进：

$$
\hat{x} \leftarrow \hat{x} + \Delta t\, F'(x)^\mathsf{T}\, \hat{x}
$$

## 3. 手写线性化（以气压梯度力为例）

动力倾向中的非线性主项为水平气压梯度力（Exner 形式）：

$$
F_u = -c_p\, \theta\, \nabla_\xi \pi, \qquad \theta = \frac{\rho\theta}{\rho}
$$

冻结背景 $\theta_b, \pi_b$，一阶线性化（扰动量带 $\delta$）：

$$
\delta F_u = -c_p \big(\delta\theta\, \nabla_\xi \pi_b + \theta_b\, \nabla_\xi
\delta\pi\big), \qquad
\delta\theta = \frac{\delta(\rho\theta) - \theta_b\,\delta\rho}{\rho_b}
$$

其伴随（转置）将伴随量反向散射：

- 经 $\delta\theta$ 项：散射到 $\delta(\rho\theta)$ 与 $\delta\rho$；
- 经 $\theta_b \nabla_\xi \delta\pi$ 项：散射到相邻网格点的 $\delta\pi$。

质量倾向 $\nabla\cdot(\rho v)$ 的线性化、科氏项（$\rho u / \rho$ 还原速度）的
线性化同理，完整推导见各算子对应的文档。

## 4. 点积检验（Dot-Product Test）

验证 TL 与 AD 严格互为转置：

$$
\frac{\big|\langle M'u, v\rangle - \langle u, M^\mathsf{T} v\rangle\big|}
{\max(|\langle M'u, v\rangle|, |\langle u, M^\mathsf{T} v\rangle|)} \lesssim 10^{-12}
$$

- 随机扰动 $u, v$ 需覆盖全状态空间（8 个字段）；
- 内积用 Kahan/长双精度补偿求和，避免大状态空间下的求和误差污染检验；
- 通过点积检验是 TL/AD 可用于 4D-Var 的**必要条件**（不充分，还需切线检验）。

## 5. 实现策略

本模块提供三层：

| 组件 | 作用 |
|---|---|
| `FiniteDifferenceTLM` | 中心差分近似 $M'$，作为手写 TL 的独立参照 |
| `HandwrittenTLM` | 手写解析线性化的单时间步传播骨架 |
| `HandwrittenADM` | 手写转置（伴随）单时间步传播骨架 |

切线检验（另需）：对任意扰动 $\delta x$，验证

$$
M(x + \varepsilon\,\delta x) - M(x) \approx M'(x)\, \varepsilon\,\delta x
$$

在 $\varepsilon \to 0$ 时线性收敛。

## 参考文献

- Errico, R. M. (1997). *What is an adjoint model?* Bull. Amer. Meteor. Soc.
  **78**, 2577–2591.
- Giering, R. & Kaminski, T. (1998). *Recipes for adjoint code construction.*
  ACM TOMS **24**, 437–474.
- Lorenc, A. C. (2003). *The potential of the ensemble Kalman filter for NWP.*
  QJRMS **129**, 3183–3203.
