# 增量 4D-Var 数据同化（P5）

> 本文描述增量 4D-Var 的成本函数、梯度计算（伴随反推）与极小化。
> 对应提示词集 P5 阶段，是 P6 切线性/伴随模型的上层消费者。

## 1. 增量 4D-Var 的代价函数

### 1.1 全量形式

4D-Var 在同化窗 $[t_0, t_N]$ 内极小化：

$$
J(x_0) = \frac12 (x_0 - x_b)^\mathsf{T} B^{-1} (x_0 - x_b)
       + \frac12 \sum_{k} \big(H_k M_k x_0 - y_k\big)^\mathsf{T} R_k^{-1}
         \big(H_k M_k x_0 - y_k\big)
$$

其中 $x_b$ 为背景（先验），$B$ 为背景误差协方差，$M_k$ 为从 $t_0$ 到观测时刻
$t_k$ 的非线性预报，$H_k$ 为观测算子，$y_k$ 为观测，$R_k$ 为观测误差协方差。

### 1.2 增量形式（Courtier et al. 1994）

为避免显式构造与求逆巨型 $B$，引入增量 $\delta x = x_0 - x_b$，并用控制变量
变换 $\delta x = B^{1/2} w$：

$$
J(w) = \frac12 w^\mathsf{T} w
     + \frac12 \sum_k \big(H_k M'_k B^{1/2} w - d_k\big)^\mathsf{T} R_k^{-1}
       \big(H_k M'_k B^{1/2} w - d_k\big)
$$

其中 $M'_k$ 为切线性模型，创新向量 $d_k = y_k - H_k M_k(x_b)$。背景项退化为
$\frac12 w^\mathsf{T} w$，天然提供了单位预条件子。

## 2. 梯度（伴随反推）

$$
\nabla_w J = w + B^{1/2\mathsf{T}} \sum_k M_k^{\prime\mathsf{T}} H_k^\mathsf{T}
R_k^{-1} \big(H_k M'_k B^{1/2} w - d_k\big)
$$

梯度计算分三步（每个时隙 $k$）：

1. **前向**：$\delta x = B^{1/2} w$，切线性传播 $\delta x_k = M'_k \delta x$；
2. **观测残差**：$r_k = H_k \delta x_k - d_k$，加权 $g_k = R_k^{-1} r_k$；
3. **伴随反推**：$\hat{x} \leftarrow H_k^\mathsf{T} g_k$，$\hat{x} \leftarrow
   M_k^{\prime\mathsf{T}} \hat{x}$，累加到 $\nabla_w J$ 的观测项。

最终 $\nabla_w J = w + B^{1/2\mathsf{T}} \hat{x}$。

## 3. 极小化（预条件共轭梯度）

由于背景项 $\frac12 w^\mathsf{T} w$ 提供单位预条件子，采用标准共轭梯度（CG）。
每步需要的 Hessian-vector 积 $\nabla^2 J \, p$ 用梯度差近似：

$$
\nabla^2 J \, p \approx \frac{\nabla J(w + \epsilon p) - \nabla J(w)}{\epsilon}
$$

（生产实现可解析计算 $\nabla^2 J$，或用 Lanczos 求解）。

## 4. 误差协方差建模

- **背景误差 $B$**：当前实现为对角（各变量独立方差 $\sigma_b^2$）。生产需
  引入空间相关（球谐谱、递归滤波、小波），参考 Bannister (2008)。
- **观测误差 $R$**：对角（各观测独立，误差标准差 $\sigma_o$）。卫星辐射率等
  需考虑通道间相关（块对角）。

## 5. 观测算子 $H$

观测算子将模式状态插值到观测位置并做变量变换（如 $\theta \to T$、湿度变量）。
本模块提供 `SimpleThetaObservationOperator` 作为接口示例（位温观测），完整实现
需接入网格度量与垂直坐标的水平/垂直插值。

## 参考文献

- Courtier, P., Thépaut, J.-N. & Hollingsworth, A. (1994). *A strategy for
  operational implementation of 4D-Var.* QJRMS **120**, 1367–1387.
- Lorenc, A. C. (2003). *The potential of the ensemble Kalman filter for NWP.*
  QJRMS **129**, 3183–3203.
- Bannister, R. N. (2008). *A review of forecast error covariance statistics.*
  QJRMS **134**, 1951–1970.
- Rawlins, F. et al. (2007). *The Met Office global four-dimensional variational
  data assimilation scheme.* QJRMS **133**, 347–362.
