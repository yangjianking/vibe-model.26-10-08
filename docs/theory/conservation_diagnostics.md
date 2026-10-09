# 守恒诊断与标准算例（P7 检验评估）

> 本文描述模式正确性检验的两大支柱：守恒诊断（离散守恒量漂移监测）与
> 标准算例（Williamson / Jablonowski / Held–Suarez 回归）。对应提示词集 P7。

## 1. 守恒诊断（Conservation Diagnostics）

### 1.1 应守恒的全球量

绝热、无摩擦、无地形的非静力可压缩 Euler 方程应精确守恒：

| 量 | 表达式 | 物理意义 |
|---|---|---|
| 总质量 | $M = \int \rho\, dV$ | 连续方程通量形式 |
| 总位温 | $\Theta = \int \rho\theta\, dV$ | 绝热无源时守恒 |
| 总能量 | $E = \int (\frac12\rho\|\mathbf v\|^2 + \rho c_v T + \rho\Phi)\, dV$ | 动能+内能+位能 |
| 总角动量 | $A = \int \rho(u\, r\cos\varphi + \Omega r^2\cos^2\varphi)\, dV$ | 轴对称性 |

离散守恒依赖散度算子与压力梯度算子的**伴随配对**（Arakawa–Lamb 1981），
保证离散总能量不产生虚假源汇。

### 1.2 全球积分的立方球实现

6 个面板的 gnomonic 等角投影积分：

$$
\int f\, dV = R^2 \sum_{p=1}^{6} \sum_{k} \int f(\xi,\eta,z_k)\, \sqrt{G}\,
d\xi\, d\eta\, \Delta z_k
$$

**关键**：`MetricPoint::sqrtG = 1/r^3` 是**单位球**（半径 1）的 Jacobian，实际
面积元需乘 $R^2$（地球半径平方）。等角坐标 $\xi,\eta\in[-1,1]$，故
$d\xi = d\eta = 2/n_{cells}$。

面积积分收敛性（Python 数值验证）：

| $n_{cells}$ | 单位球面积 | 相对误差 |
|---|---|---|
| 8 | 12.61449335 | 3.829e-3 |
| 16 | 12.57839964 | 9.572e-4 |
| 48 | 12.56770709 | 1.064e-4 |
| 96 | 12.56670473 | 2.659e-5 |

面积积分随分辨率一阶收敛到 $4\pi$，与理论一致。

### 1.3 补偿求和

全球积分涉及 $\sim 6 n^2 n_{lev}$ 个浮点累加，朴素求和会积累舍入误差
（相对误差 $\sim \sqrt{N}\,\epsilon$）。采用**长双精度（80-bit extended）累加器**
避免求和误差污染守恒诊断（提示词 P3 禁朴素求和）。

## 2. 标准算例（Test Cases）

### 2.1 Williamson TC2（稳态定常平流）

余弦钟型标量场被恒定角速度的定常风场平流（12 天绕地球一周），**有解析解**，
用于误差范数与收敛阶检验。

- 速度场：$u = u_0(\cos\varphi\cos\alpha + \sin\varphi\cos\lambda\sin\alpha)$，
  $v = -u_0\sin\lambda\sin\alpha$，$u_0 = 2\pi a / 12\text{天}$。
- 钟型场：$h(r) = \frac12[1 + \cos(\pi r / R)]$，$r < R = \pi/3$。
- 解析解：场以角速度 $\Omega = u_0/a$ 沿 $\alpha$ 方向刚体平移。

### 2.2 Williamson TC5（地形罗斯贝波）

恒定纬向流跨越孤立高斯山，产生定常罗斯贝波。地形高度
$h_s = h_0(1 - r/R)$，$h_0 = 2000$ m，$R = \pi/9$。检验地形处理与平衡性。

### 2.3 Jablonowski–Williamson 斜压波

解析的静力平衡斜压波初值（JW2006），含纬向喷流与扰动，检验三维非静力
动力核心的斜压响应。

### 2.4 Held–Suarez 气候态强迫

理想化气候态的**强迫**（非初值算例）：温度弛豫到解析平衡态 + 边界层瑞利摩擦。

- 平衡温度（HS94 Eq. 3）：
  $$
  T_{eq}(\varphi,\sigma) = \max\Big\{200,\;
  \big[315 - 60\sin^2\varphi - 10\cos^2\varphi\,\tfrac{\sigma-\sigma_b}{1-\sigma_b}\big]
  \,\sigma^{\kappa}\Big\},\quad \sigma_b = 0.7
  $$
- 瑞利摩擦（HS94 Eq. 5）：$k_v = \frac{1}{1\text{天}}\frac{\sigma-\sigma_b}{1-\sigma_b}$
  for $\sigma > \sigma_b$。

> 注：本项目用高度坐标，$\sigma = p/p_s \approx 1 - z/z_{top}$（低层 $\sigma\to1$）。

## 3. 检验流程

1. **守恒漂移**：长期积分中监测 $M, \Theta, E, A$ 相对初值的漂移，应达
   机器精度（离散守恒）或缓慢收敛。
2. **TC2 收敛阶**：不同分辨率下计算 L1/L2/Linf 误差范数，验证平流方案收敛阶。
3. **TC5 平衡**：长时间积分后罗斯贝波解应保持定常。
4. **HS94 气候态**：多年积分后达到统计平衡态（纬向平均温度/风场）。

## 参考文献

- Arakawa & Lamb (1981), Mon. Wea. Rev. **109**, 18–36.
- Williamson et al. (1992), J. Comput. Phys. **102**, 211–224.
- Held & Suarez (1994), Bull. Amer. Meteor. Soc. **75**, 1825–1830.
- Jablonowski & Williamson (2006), QJRMS **132**, 2943–2975.
- Thuburn (2008), J. Comput. Phys. **227**, 3715–3730.
