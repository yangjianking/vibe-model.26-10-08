# 未完成工作的数学推导问题集（Mathematical Problems to Be Derived）

> 本文档将 `TODO.md` 中所有**需要数学推导**的未完成工作，按严格的数学问题形式
> （问题陈述 → 已知/假设 → 待证结论 → 推导路径 → 期望结果）逐一书写。
> 与 `TODO.md` 的工程清单互补：TODO 记"做什么"，本文档记"怎么从数学上推出来"。
> 对应提示词集 P0 第 6 条「连续形式 → 离散化推导 → 截断误差阶数 → 稳定性约束」。
>
> 最后更新：2026-10-09

---

## 问题 1 · 半隐式 Helmholtz 方程的完整消元推导

**对应**：TODO #6（Helmholtz 多重网格预条件）、`docs/theory/helmholtz.md` §2（当前为概要）。

**已知**：非静力可压缩 Euler 方程的线性化快波系统（在参考态 $(\bar\rho,\bar\theta,\bar\pi)$ 附近）：

$$
\frac{\partial \delta\mathbf v}{\partial t}
= -c_p\bar\theta\,\nabla(\delta\pi') - g\frac{\delta\rho}{\bar\rho}\hat{\mathbf k} + \cdots
$$

$$
\frac{\partial \delta\pi'}{\partial t}
= -\frac{c_s^2}{c_p\bar\theta}\,\nabla\cdot(\bar\rho\,\delta\mathbf v) + \cdots
$$

其中 $\delta\pi'$ 为 Exner 压力扰动，$c_s^2 = \gamma R_d\bar T$（参考态声速平方），
$\hat{\mathbf k}$ 为垂直单位向量，$\cdots$ 为慢过程与非线性项（显式处理）。

**待证**：对快波项做 Crank–Nicolson 隐式平均（off-centering 参数
$\alpha\in[0.5,0.6]$，权重 $\beta=\alpha\Delta t$），消去 $\delta\mathbf v$ 后，
得到关于 $\delta\pi'$ 的三维椭圆方程：

$$
\Big(I - \beta^2\Delta t^2\, \nabla\cdot\!\big(\bar\rho\, c_s^2\, \nabla\big)\Big)\delta\pi'
= R
$$

其中右端项 $R$ 由显式中间速度场 $\mathbf v^*$ 的散度给出：

$$
R = -\beta\Delta t\,\nabla\cdot(\bar\rho\,\mathbf v^*)
$$

**推导路径**：
1. 将动量方程（仅快波线性项）在 $[t_n,t_{n+1}]$ 上用 Crank–Nicolson 离散，
   解出 $\delta\mathbf v^{n+1} = \delta\mathbf v^* - \beta\Delta t\, c_p\bar\theta\nabla(\delta\pi')$，
   其中 $\delta\mathbf v^*$ 为显式预报的中间速度；
2. 将 $\delta\mathbf v^{n+1}$ 代入连续方程的隐式部分，得到关于 $\delta\pi'$ 的
   单一椭圆方程；
3. 明确 $\nabla\cdot(\bar\rho\,\nabla)$ 是**自伴**（对称正定）算子，$I - \beta^2\Delta t^2$
   乘该算子是**对称正定**的，故 BiCGStab/GMRES 收敛有理论保证；
4. 说明为何在深大气下 $c_s$ 随高度变化时算子仍保持自伴（权重 $\bar\rho c_s^2$ 在
   散度内侧）。

**期望结果**：一条从线性化方程到 Helmholtz 方程的**逐步代数推导**，明确每一步
丢弃的项及其量级，并给出算子对称正定的证明。

**文献锚点**：Wood et al. (2014, QJRMS 140)；Staniforth & Côté (1991, MWR 119)；
Tanguay et al. (1990, MWR 118)。

---

## 问题 2 · 半隐式格式的完整 von Neumann 稳定性分析

**对应**：`docs/numerics/semi_implicit.md` §4（当前仅给放大因子结论）。

**已知**：线性化快波系统（问题 1 的连续形式），对单波解
$q \propto e^{i(\mathbf k\cdot\mathbf x - \omega t)}$ 做模态分析。

**待证**：
1. **Crank–Nicolson 半隐式对重力波/声波无条件稳定**：放大因子
   $A = (1+i\alpha\omega\Delta t)/(1-i\alpha\omega\Delta t)$ 对任意实 $\omega$ 满足
   $|A|=1$（中性稳定，无耗散）；
2. **off-centering 参数 $\alpha > 0.5$ 引入数值阻尼**，阻尼率 $\rho$ 的显式表达式为：

$$
\rho = \frac{1}{1 + (2\alpha-1)^2\omega^2\Delta t^2}
\qquad\text{（对 $2\Delta t$ 计算模态）}
$$

并证明 $\alpha=0.5$（纯 Crank–Nicolson）时 $2\Delta t$ 模态无阻尼（弱不稳定），
故需 $\alpha\in(0.5,0.6]$；
3. **平流项**对显式部分给出 CFL 约束 $\Delta t_{adv} \le C_{\max}\cdot(\Delta x/u_{\max})$，
   证明半隐式消除声波 CFL 后，$\Delta t$ 仅受平流精度约束；
4. 三维耦合情形下，给出**垂直传播重力波**的色散关系与最大可稳定 $\Delta t$ 估计。

**推导路径**：对每个波模代入离散方程，求放大因子矩阵的特征值，讨论谱半径。

**期望结果**：完整的放大因子推导（非概要），含 $2\Delta t$ 模态阻尼率公式与
$\alpha$ 取值下限的证明。

**文献锚点**：Robert (1981, Atmos.-Ocean 19)；Staniforth & Côté (1991)；
Durran (2010, 第 2 版第 3 章 von Neumann 分析)。

---

## 问题 3 · 离散总能量守恒的代数证明（Arakawa–Lamb 配对）

**对应**：`docs/theory/dynamics_equations.md` §4.2（当前为"概要，后续补全"）、
`docs/numerics/pressure_gradient.md` §4。

**已知**：离散的动量方程（协变分量）与热力学方程（位温密度形式），以及散度算子
与压力梯度算子的离散模板（C 网格交错，立方球度量）。

**待证**：若散度算子 $\nabla_d\cdot$ 与压力梯度算子 $\nabla_d$ 满足**离散伴随关系**：

$$
\sum_{i,j} \mathbf v_{i,j}\cdot(\nabla_d \pi)_{i,j}\,\sqrt G\,\Delta\xi\Delta\eta
= -\sum_{i,j} \pi_{i,j}\,(\nabla_d\cdot\mathbf v)_{i,j}\,\sqrt G\,\Delta\xi\Delta\eta
$$

（即两者互为离散转置/负伴随），则离散总能量

$$
E = \sum \Big(\tfrac12 \rho|\mathbf v|^2 + \rho c_v T + \rho gz\Big)\sqrt G\,\Delta\xi\Delta\eta\,\Delta z
$$

满足 $\dfrac{dE}{dt}=0$（无源汇时）。

**推导路径**：
1. 动量方程点乘 $\mathbf v$（协变点乘需经度量张量 $g_{ij}$），对全球求和；
2. 热力学方程乘 $c_p\theta$（或对能量形式用 $c_v T$），对全球求和；
3. 两式相加，利用散度/梯度伴随配对使气压梯度做功项与散度项精确抵消；
4. 证明**跨面板通量项**经面板连续性与 halo 交换后望远镜式求和为零
   （flux pairing 守恒，已在平流算子中实测机器精度）；
5. 给出离散角动量 $A$ 与位温 $\Theta$ 的守恒条件（轴对称与无源）。

**期望结果**：逐步代数证明，明确每一步依赖的离散结构（交错位置、度量、通量配对），
并指出哪些项会破坏守恒（如耗散加热是否回补，Harris et al. 2021 §8.5）。

**文献锚点**：Arakawa & Lamb (1981, MWR 109)；Arakawa (1966, JCP 1)；Thuburn (2008, JCP 227)。

---

## 问题 4 · 几何多重网格的局部 Fourier 分析（LFA）

**对应**：TODO #6（多重网格预条件）、`docs/theory/helmholtz.md` §4（仅预留）。

**已知**：三维椭圆算子 $A = I - \beta^2\Delta t^2\, c_s^2\nabla^2$（立方球各向异性
度量 $g^{ij}$，垂直地形跟随坐标强耦合）。

**待证**：
1. 离散算子 $A$ 的**谱半径**上下界估计：
   $\lambda_{\min}(A) \ge 1$，$\lambda_{\max}(A) \le 1 + \beta^2\Delta t^2 c_s^2\,\lambda_{\max}(\nabla^2)$，
   并给出 $\lambda_{\max}(\nabla^2)$ 在立方球网格上的显式估计
   （$\sim O(1/\Delta x^2)$ 量级，各向异性系数 $\max g^{ij}$）；
2. 对**点 Jacobi / Gauss-Seidel 光滑器**做 LFA，给出光滑因子 $\mu$ 与各向异性比
   $\epsilon = g^{11}/g^{22}$ 的显式关系，证明当 $\epsilon \gg 1$（面板边缘各向异性
   强烈）时点光滑退化，需**线光滑/块光滑**；
3. 给出两网格法的收敛因子 $\rho_{2h}$ 估计，说明为何垂直方向需**线松弛/ADI**
   （垂直强耦合），水平周期方向可用 FFT。

**推导路径**：LFA 标准框架——将光滑器作用于单个 Fourier 模
$e^{i\mathbf k\cdot\mathbf x}$，求其放大因子，取高/低频空间的最大值。

**期望结果**：给出光滑因子随各向异性比的变化曲线（解析或数值），论证多重网格
预条件在立方球上的可行性，并确定光滑器选型。

**文献锚点**：Buckeridge & Scheichl (2010, QJRMS 136)；Trottenberg et al. (2001, *Multigrid*)；
Dedner et al. (2016, GMD 9)。

---

## 问题 5 · 背景误差协方差 B 的建模与控制变量变换（CVT）

**对应**：TODO #8（背景误差相关模型）、`docs/theory/four_dvar.md` §4（当前为对角占位）。

**已知**：背景误差协方差 $B \in \mathbb R^{n\times n}$（对称正定，$n$ 为控制空间维数，
可达 $10^8$ 量级），无法显式构造/存储。

**待证**：
1. **控制变量变换** $\delta x = B^{1/2} w$ 将背景项
   $\tfrac12\delta x^\mathsf T B^{-1}\delta x$ 化为 $\tfrac12 w^\mathsf T w$（单位
   预条件），需证明 $B^{1/2}$（对称平方根）的存在性与作用方式（矩阵-向量积）；
2. **递归滤波**算子 $L$ 满足 $B \approx L L^\mathsf T$，证明一阶自回归（AR(1)）
   递归滤波的协方差结构为各向同性指数相关
   $\rho(r) = e^{-r/L_c}$（$L_c$ 为相关长度），并给出其谱表示；
3. **球谐谱滤波**：水平均匀各向同性假设下，$B$ 的谱分解
   $B = \sum_\ell b_\ell\, P_\ell$（$P_\ell$ 为 $\ell$ 阶谱投影），给出方差谱
   $b_\ell$ 与垂直递归相关的张量积结构；
4. **平衡算子**：$\delta x = K\,\delta x_u + \delta x_{bal}$ 中，平衡算子 $K$
   （地转/统计回归）的定义与最小二乘标定，证明其将变量间的平衡约束编码进 $B$。

**推导路径**：从 $B$ 的假设结构（同质各向同性水平相关 × 垂直相关 × 变量平衡耦合）
出发，逐步构造可作用的 $B^{1/2}$，并对每步给出作用一次的算法复杂度。

**期望结果**：$B^{1/2}$ 的三种构造（递归滤波/谱/小波）的完整数学推导与复杂度
对比，明确本项目首选（建议递归滤波，纯本地可实现）。

**文献锚点**：Bannister (2008, QJRMS 134)；Parrish & Derber (1992, MWR 120)；
Fisher (2003, ECMWF TM)；Derber & Bouttier (1999, Tellus 51A)。

---

## 问题 6 · 切线性模型的切线检验与伴随的严格转置性

**对应**：TODO #7（手写 TL/AD 切线检验）、`docs/theory/tlm_ad.md` §4-5。

**已知**：预报模型 $M: x_0 \mapsto x_N$ 及其手写切线性 $M'$、手写伴随 $M^{\prime\mathsf T}$。

**待证**：
1. **切线检验**（Taylor 余量）：对任意扰动方向 $d$，比值

$$
\phi(\varepsilon) = \frac{\|M(x+\varepsilon d) - M(x)\|}{\|\varepsilon M'(x)d\|}
$$

满足 $\phi(\varepsilon) \to 1$ 当 $\varepsilon\to 0$，且误差
$|1 - \phi(\varepsilon)| = O(\varepsilon)$（一阶线性收敛），证明 $M'$ 是 $M$ 的
真切线性（Gateaux 导数）；
2. **点积检验**（转置性，已实测通过 $2.65\times10^{-16}$）：证明
   $\langle M'u, v\rangle = \langle u, M^{\prime\mathsf T}v\rangle$ 对任意 $u,v$ 成立
   当且仅当 ADM 是 TLM 的严格代数转置（"先离散后线性化"原则）；
3. **非光滑点处理**：对 $q = \max(0,\cdot)$、$\operatorname{sign}(\cdot)$、条件跳转
   处的线性化，证明子梯度选择（或正则化光滑替代）对梯度精度的误差界。

**推导路径**：对每个动力算子（平流通量、压力梯度、半隐式求解器、重映射），写
连续方程 → 扰动方程（TL）→ 转置离散（AD）三段式推导，再整体做切线/点积检验。

**期望结果**：每个算子的 TL/AD 三段式推导 + 切线检验收敛曲线（$O(\varepsilon)$）+ 
点积检验机器精度报告。

**文献锚点**：Giering & Kaminski (1998, ACM TOMS 24)；Errico (1997, BAMS 78)；
Talagrand & Courtier (1987, QJRMS 113)；Janisková et al. (1999, Tellus 51A)。

---

## 问题 7 · 4D-Var 代价函数梯度的伴随表达式与增量式等价性

**对应**：TODO #9（多时隙完整 4D-Var）、`docs/theory/four_dvar.md` §2（梯度已给，缺等价性论证）。

**已知**：增量 4D-Var 代价函数（问题 5 的 CVT 后）：

$$
J(w) = \tfrac12 w^\mathsf T w + \tfrac12\sum_k \big(H_k M'_k B^{1/2}w - d_k\big)^\mathsf T
R_k^{-1}\big(H_k M'_k B^{1/2}w - d_k\big)
$$

**待证**：
1. 梯度表达式
   $\nabla_w J = w + B^{1/2\mathsf T}\sum_k M_k^{\prime\mathsf T} H_k^\mathsf T R_k^{-1}(H_k M'_k B^{1/2}w - d_k)$
   的**逐项推导**（对每个时隙的观测项，用链式法则 + 伴随反推）；
2. **增量式与扩展 Kalman 滤波（EKF）在窗末的等价性**（Courtier et al. 1994 §3）：
   证明增量式 4D-Var 的分析增量 $\delta x_a$ 与 EKF 在窗末的更新等价（在线性、
   高斯、完美模式假设下）；
3. **多增量（multi-incremental）外循环**的收敛性：外循环用高分辨率非线性轨迹
   更新新息，内循环低分辨率 TLM/ADM 极小化，证明该分裂的合理性；
4. Hessian 的预条件：证明背景项 $\tfrac12 w^\mathsf T w$ 提供**单位预条件子**，
   使内循环 CG 收敛仅依赖观测项的谱条件数。

**期望结果**：梯度逐项推导 + 增量式与 EKF 等价性论证 + 内外循环收敛性说明。

**文献锚点**：Courtier et al. (1994, QJRMS 120)；Veersé & Thépaut (1998, QJRMS 124)；
Rabier et al. (2000, QJRMS 126)。

---

## 问题 8 · 切线性物理的线性化（降水/边界层/扩散）

**对应**：`docs/theory/tlm_ad.md` §3（仅给气压梯度例）、TODO 中物理参数化（#5）。

**已知**：物理参数化方案（大尺度降水、边界层拖曳、垂直扩散）是非线性算子
$P: q \mapsto P(q)$（含开关/阈值）。

**待证**：对每一类物理过程，给出其**线性化版本** $P'(q)$（用于同化，参考
Janisková & Lopez 2013 的"同化用线性物理"包）：

1. **大尺度降水**：凝结率 $C(q_v, T)$ 的线性化 $\delta C = \partial_{q_v}C\,\delta q_v +
   \partial_T C\,\delta T$，给出偏导数的显式形式（Clausius–Clapeyron 关系）；
2. **边界层拖曳**：$F = -c_d|\mathbf v|\mathbf v$ 的线性化
   $\delta F = -c_d(|\mathbf v|\,\delta\mathbf v + \mathbf v\,\delta|\mathbf v|)$；
3. **垂直扩散**：$\partial_z(K\partial_z q)$ 的线性化（$K$ 冻结或线性化）；
4. 对开关/阈值（降水开关、云开关）给出子梯度或正则化光滑替代，并分析其对
   梯度精度的误差界。

**期望结果**：每个物理算子"连续 → 扰动方程（TL）→ 转置（AD）"三段式，及开关处理
的误差分析。

**文献锚点**：Janisková et al. (1999, Tellus 51A)；Mahfouf (1999, QJRMS 125)；
Janisková & Lopez (2013, ECMWF TM)。

---

## 问题 9 · 嵌套边界的守恒性与反射分析

**对应**：提示词 P1（嵌套，one-way/two-way）、TODO 中跨面板 halo（#2, #10）。

**已知**：粗细网格比 1:3，粗→细用 Davies (1976) 松弛 + 海绵层，细→粗用通量形式
保守替换（flux-form replacement）。

**待证**：
1. **Davies 松弛**：松弛系数 $\gamma(r)$ 随距离 $r$ 的函数形式，及反射系数
   $|R(\omega)|$ 的分析——证明对目标波频段的反射振幅 $< 2\%$；
2. **海绵层**：阻尼系数 $\mu(r)$ 的衰减函数，及对进入海绵层的波的能量吸收率；
3. **通量形式保守替换**：细网格积分量 $\int f\,dV$ 与粗网格积分量在嵌套界面
   的守恒性证明（质量/能量守恒误差 $< 1\times10^{-10}$ 相对）；
4. **跨面板 halo 交换**：向量场跨面板的旋转矩阵 $R$（问题见 `cubed_sphere_metric.md`
   §6）与通量一致性，证明面板边界通量的望远镜式求和为零（全球守恒）。

**期望结果**：松弛系数/海绵层函数形式的完整推导 + 反射系数分析 + 守恒性证明。

**文献锚点**：Davies (1976, QJRMS 102)；Harris & Lin (2013, MWR 141)；
Schmidt (1977, Beitr. Phys. Atmos. 50)。

---

## 问题 10 · Schär 地形跟随坐标的度量与坡度限制

**对应**：`docs/theory/vertical_coordinate.md`（度量已给，坡度限制待补）、
TODO 中地形处理。

**已知**：Schär (2002) 平滑地形跟随坐标 $z(\eta) = z_{top}\,\eta + h_s\,b(\eta)$，
$b(\eta) = \dfrac{\sinh\!\big((1-\eta)/s\big)}{\sinh(1/s)}$（$s$ 为衰减参数）。

**待证**：
1. 垂直 Jacobian $J = \partial z/\partial\eta$ 与度量项 $\partial z/\partial x$、
   $\partial z/\partial y$ 的完整推导，及水平导数算子在坐标变换下的形式
   （链式法则，$\nabla_x = \nabla_\eta - \frac{1}{J}\nabla_\eta z\,\partial_\eta$）；
2. **单调性条件**：证明 $J > 0$ 的充分条件（坐标面不交叉），并给出地形高度
   $h_s$ 与衰减参数 $s$ 的允许范围；
3. **Janjić (1984) 坡度限制**：给出地形坡度 $\gamma = \arctan|\nabla h_s|$ 的
   上限判据（避免 $J$ 过小导致的坐标退化），及其与格点分辨率的显式关系；
4. 证明该坐标在静力平衡下的压力梯度力离散不产生虚假环流（地形截断误差）。

**期望结果**：坐标变换的完整度量推导 + 单调性/坡度限制的显式判据。

**文献锚点**：Schär et al. (2002, MWR 130)；Janjić (1984)；Simmons & Burridge (1981)。

---

## 问题 11 · 平流方案截断误差阶数（modified equation 分析）

**对应**：`docs/numerics/advection.md`（PPM 已给，缺阶数推导）、TODO #11（收敛阶）。

**已知**：三种平流方案——中心二阶、三阶上游（Wicker–Skamarock 2002）、PPM
（Colella & Woodward 1984）。

**待证**：对一维线性平流方程 $\partial_t q + u\partial_x q = 0$，用 modified equation
分析给出各方案的**有效截断误差**：

1. 中心二阶：主导误差项 $O(\Delta x^2)$ 为**色散**（相位误差），无耗散；
2. 三阶上游：主导误差项 $O(\Delta x^3)$ 为**耗散**（数值扩散系数显式给出）；
3. PPM（含单调限制器）：理想 $O(\Delta x^4)$ 插值，但单调限制器在极值处降阶，
   解释实测收敛阶 2.35（介于 2 与 3 之间，因限制器退化）；
4. 给出 SSP-RK3 时间推进下各方案的**综合收敛阶**（空间+时间误差耦合），与
   `python/cubedsphere/convergence.py` 实测（2.00 / 2.98 / 2.35）对照。

**期望结果**：各方案的 modified equation 主导项 + 与实测收敛阶的定量对照。

**文献锚点**：Colella & Woodward (1984, JCP 54)；Wicker & Skamarock (2002, MWR 130)；
LeVeque (2002, 第 2 版第 8-9 章)。

---

## 问题 12 · 标准算例解析解/平衡解的推导

**对应**：TODO #11–14（TC2/TC5/HS94 长时间积分）、`docs/theory/conservation_diagnostics.md` §2。

**待证**：
1. **Williamson TC2 解析解**：证明余弦钟型场 $h$ 在恒定角速度风场下做**刚体平移**
   （解析解 $h(\mathbf x, t) = h_0(\mathbf x - \Omega t)$），给出其任意时刻的解析
   表达式（用于误差范数）；
2. **Williamson TC5 定常解**：证明恒定纬向流跨越高斯山的稳态解是**定常罗斯贝波**
   （位涡守恒 $\frac{D}{Dt}(\zeta + f)/h = 0$），给出定常波形的解析或参考形式；
3. **Held–Suarez 平衡温度**：推导 $T_{eq}(\varphi,\sigma)$ 的构造（HS94 Eq. 3），
   证明 $\sigma\to 1$（低层）时 $T_{eq}$ 趋于热带 ~315 K、极区 ~200 K 的合理性；
4. **Jablonowski–Williamson 初值**：给出解析的静力平衡斜压波初值（纬向喷流 +
   扰动），证明其满足静力平衡与热成风关系。

**期望结果**：各算例的解析解/平衡解/初值构造的完整推导，作为检验的"真值"。

**文献锚点**：Williamson et al. (1992, JCP 102)；Held & Suarez (1994, BAMS 75)；
Jablonowski & Williamson (2006, QJRMS 132)。

---

## 问题 13 · 全球积分与补偿求和的误差界

**对应**：`docs/theory/conservation_diagnostics.md` §1.3（已给，可补误差界证明）、
TODO #14（长期漂移）。

**已知**：全球积分 $S = \sum_{m=1}^N a_m$（$N \sim 6 n_{cells}^2 n_{lev}$ 项，
浮点累加）。

**待证**：
1. **朴素求和的误差界**：$\big|S - \hat S\big| \le \gamma_{N-1}\sum_m |a_m|$，
   其中 $\gamma_n = n\varepsilon/(1-n\varepsilon)$（$\varepsilon$ 机器精度），
   说明误差随 $N$ 线性增长；
2. **Kahan 补偿求和**：误差界为 $O(\varepsilon + N\varepsilon^2)$，与 $N$ 无关，
   证明补偿项 $c$ 捕获了低位舍入；
3. **长双精度（80-bit）累加器**：$\varepsilon_{80} \approx 5.4\times10^{-20}$，
   对 $N \sim 10^9$ 项误差仍可忽略；
4. 给出守恒诊断相对误差 $< 1\times10^{-10}$ 所需的累加精度选择。

**期望结果**：三类求和的误差界证明与选型建议（已选长双精度，补证明）。

**文献锚点**：Higham (2002, *Accuracy and Stability*, 第 4 章)。

---

## 问题 14 · 混合精度舍入误差的误差预算

**对应**：`docs/numerics/mixed_precision.md`（已给定性论证，可补定量证明）、TODO 中
混合精度长积分漂移。

**已知**：主体 fp32（$\varepsilon_{32}=1.19\times10^{-7}$），关键量 fp64。

**待证**：
1. 单步平流/物理倾向的相对舍入误差 $O(\varepsilon_{32})$，$N$ 步无偏累积后
   $O(\sqrt N\,\varepsilon_{32})$（随机游走模型），对 10 天（$N\sim 10^4$）为
   $\sim 10^{-5}$；
2. 证明 $\sim 10^{-5}$ 远小于离散误差（$10^{-3}\sim10^{-4}$）与物理不确定性
   （$10^{-2}\sim10^{-1}$），故 fp32 主体**不主导**误差预算；
3. **哪些量必须 fp64**：给出网格几何（$g_{ij},\sqrt G$）、椭圆求解器残差、全局
   守恒诊断、状态累积量降为 fp32 会破坏守恒/位可重复的**定量论证**（守恒量对
   舍入的敏感性）；
4. 位可重复条件：证明同精度同分区两次运行的 bit-identical 需要确定性归约
   （固定顺序、无 atomicAdd），给出确定性归约树的误差性质。

**期望结果**：误差预算的定量证明（非定性），含各降精度点的误差传播分析。

**文献锚点**：Váňa et al. (2017, MWR 145)；Düben & Palmer (2014, MWR 142)；
Higham (2002)；Zhang et al. (2024, GMD 17)。

---

## 汇总表：问题 → 对应 TODO → 依赖 → 文献

| 问题 | 对应 TODO 项 | 主要文献（ima 可检索） | 推导难度 |
|---|---|---|---|
| 1 Helmholtz 消元 | #6 | Wood 2014、Staniforth & Côté 1991 | ★★★ |
| 2 von Neumann 稳定性 | #6 | Robert 1981、Durran 2010 | ★★★ |
| 3 离散能量守恒 | #10 | Arakawa & Lamb 1981、Thuburn 2008 | ★★★★ |
| 4 多重网格 LFA | #6 | Buckeridge 2010、Trottenberg 2001 | ★★★★ |
| 5 B 建模与 CVT | #8 | Bannister 2008、Fisher 2003 | ★★★ |
| 6 TL/AD 切线+转置 | #7 | Giering 1998、Errico 1997 | ★★★ |
| 7 4D-Var 梯度+等价性 | #9 | Courtier 1994、Veersé 1998 | ★★★ |
| 8 切线性物理 | #5 | Janisková 1999、Mahfouf 1999 | ★★★ |
| 9 嵌套守恒+反射 | #2, #10 | Davies 1976、Harris & Lin 2013 | ★★★ |
| 10 Schär 度量+坡度 | 地形 | Schär 2002、Janjić 1984 | ★★ |
| 11 平流阶数 | #11 | Colella & Woodward 1984、LeVeque 2002 | ★★ |
| 12 算例解析解 | #11–14 | Williamson 1992、Held & Suarez 1994 | ★★ |
| 13 补偿求和误差界 | #14 | Higham 2002 | ★★ |
| 14 混合精度误差预算 | P3 | Váňa 2017、Higham 2002 | ★★ |

---

## 说明

- 上述问题**全部**为"接口已就绪、数学推导待补全"的实质内容，与 TODO.md 的
  工程清单一一对应。
- 其中**问题 2、3、4、5、6、8** 是提示词集 P0 第 6 条明令"禁止跳步"的核心推导，
  优先级最高；**问题 1、7、9、10、11、12、13、14** 是各阶段验收标准的数学依据。
- 所有文献锚点均可在 ima 知识库检索到原文（详见 `KNOWLEDGE_GAP_ANALYSIS.md`）。
