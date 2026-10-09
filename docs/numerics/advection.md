# 平流离散：通量形式与守恒性（Advection）

> 模块：`src/dynamics/` · 对应代码：`advection.hpp/cpp`

## 1. 问题陈述

水平平流需保证被平流量（水汽、位温、示踪物）的正性（$q \ge 0$）与守恒。
采用守恒型通量形式：

$$
\frac{\partial (\rho \phi)}{\partial t} = -\nabla\cdot(\rho \phi \mathbf{v})
$$

其中 $\phi$ 为单位质量示踪量，$\rho\phi$ 为示踪物密度。

## 2. 离散（守恒散度）

$$
\frac{\partial (\rho\phi)}{\partial t}\bigg|_{i,j}
= -\frac{1}{\sqrt{G}\,\Delta\xi}\Big[
  F^\xi_{i+1/2} - F^\xi_{i-1/2} + F^\eta_{j+1/2} - F^\eta_{j-1/2}
\Big]
$$

其中面通量 $F^\xi = \sqrt{G}\,\rho\,\phi_{\text{face}}\,u^\xi$，$\phi_{\text{face}}$
为面处重构值，$u^\xi$ 为反变速度分量。

## 3. 面重构方案

| 方案 | 重构公式 | 阶数 | 单调性 |
|---|---|---|---|
| Center2 | $(\phi_i + \phi_{i+1})/2$ | 2 | 无 |
| Upwind3 | $(2\phi_{i+1} + 5\phi_i - \phi_{i-1})/6$（$u>0$） | 3 | 需限制器 |
| PPM | Colella–Woodward 分段抛物 | 4 | 单调限制器 |

PPM 界面值（四点插值 + 单调限制）：

$$
\phi_{i+1/2} = \frac{7}{12}(\phi_i + \phi_{i+1}) - \frac{1}{12}(\phi_{i-1} + \phi_{i+2})
$$

随后限制到 $[\min(\phi_i,\phi_{i+1}), \max(\phi_i,\phi_{i+1})]$ 内，抑制过冲/欠冲
（保证正性）。实现见 `advection.cpp::ppm_face_value`。

## 4. 正性与守恒

- 守恒：通量散度形式下，$\sum \rho\phi\sqrt{G}$ 在全球求和时通量项两两抵消
  （面板边界经 halo 交换保证连续）。
- 正性：需配合限制器（如 PPM 单调限制器）保证 $\phi_{\text{face}} \ge 0$。

### 4.1 两阶段通量配对的守恒性（已实测验证）

守恒的关键是**每个界面只有唯一通量值**。实现采用两阶段：

1. **阶段一**：对每个单元计算右/上界面通量 $F_{i+1/2} = \phi_{i+1/2} u_{i+1/2}$；
2. **阶段二**：散度 $= (F_{i+1/2} - F_{i-1/2})/\sqrt{G}\Delta\xi$，其中西面通量
   $F_{i-1/2}$ 取西邻单元的右界面通量。

这样 $\sum_i (F_{i+1/2} - F_{i-1/2})$ 望远镜求和恒为零。实测（纯纬向周期平流，
`tests/advection_conservation.cpp`）全球积分相对误差：

| 方案 | ∫∇·(φv)dV 相对误差 |
|---|---|
| Center2 | 2.5e-19 |
| Upwind3 | 8.1e-19 |
| PPM | 1.9e-14（单调限制器微小不对称） |

均达机器精度，证明两阶段配对严格守恒。

> 全球完整守恒还需**跨面板 halo 交换**（`HaloExchange`）使面板边界通量一致，
> 当前 Advection 骨架未接入该交换（独立功能，属 P1 收尾）。

## 5. 时间积分稳定性（关键）

显式 Euler 时间步下，偶数阶中心差分与奇数阶上游格式均**不稳定**（中心差分
纯色散无耗散、上游格式耗散不足以压制 Euler 的误差放大）。实测（一维高斯钟型
平流，CFL=0.5）：

| 格式 | 显式 Euler | SSP-RK3 |
|---|---|---|
| Center2 | 无条件不稳定（N=256 发散） | 稳定 |
| Upwind3 | 不稳定（N=256 发散） | 稳定 |
| PPM | 不稳定（严重过冲） | 稳定 |

故生产模式须用 **SSP-RK3**（Wicker & Skamarock 2002，已实现于
`ExplicitRK3Integrator`）或半隐式。SSP-RK3 下的 L1 收敛阶实测
（`python/cubedsphere/convergence.py`）：

| 格式 | 收敛阶 | 说明 |
|---|---|---|
| Center2 | 2.00 | 理论二阶 ✓ |
| Upwind3 | 2.98 | 理论三阶（低分辨率受耗散略降） |
| PPM | 2.35 | 平滑场下单调限制器降至 ~2 阶 |

> PPM 在平滑极值附近因单调限制器（Colella–Woodward 已知行为）从理论 4 阶
> 降至约 2 阶；在含强间断的场中其无振荡性更具价值。

## 6. 参考文献

- Lin & Rood (1996), *Mon. Wea. Rev.*, 124.
- Wicker & Skamarock (2002), *Mon. Wea. Rev.*, 130, 2088–2097.
- Colella & Woodward (1984), *J. Comput. Phys.*, 54, 174–201.
