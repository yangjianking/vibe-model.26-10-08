# 科氏力与曲率项（Coriolis & Curvature）

> 模块：`src/dynamics/` · 对应代码：`coriolis.hpp/cpp`

## 1. 问题陈述

球面上动量方程的科氏力与曲率项需在立方球协变分量表示下正确离散。采用
**协变分量**存储风速的关键优势（Putman & Lin 2007 §3）：相对涡度可直接由
协变分量计算，避免度量 Christoffel 符号的显式处理。

## 2. 协变形式的科氏力与曲率项

协变动量方程中，科氏力与曲率项（含动能梯度）为：

$$
\frac{\partial u_\xi}{\partial t} \leftarrow (f + \zeta)\, u_\eta
- \frac12 \frac{\partial}{\partial\xi}\left(g^{ij}u_i u_j\right)
$$

$$
\frac{\partial u_\eta}{\partial t} \leftarrow -(f + \zeta)\, u_\xi
- \frac12 \frac{\partial}{\partial\eta}\left(g^{ij}u_i u_j\right)
$$

其中：
- $u_\xi, u_\eta$：协变风速分量
- $f = 2\Omega\sin\phi$：科氏参数（$\Omega = 7.2921\times10^{-5}$ rad/s）
- $\zeta$：相对涡度
- $g^{ij}u_i u_j$：动能（用逆度量张量）

## 3. 相对涡度（协变分量直接计算）

$$
\zeta = \frac{1}{\sqrt{G}}\left(
\frac{\partial u_\eta}{\partial\xi} -
\frac{\partial u_\xi}{\partial\eta}
\right)
$$

采用协变分量后，涡度算子不含显式 Christoffel 项（Putman & Lin 2007），
对应 `coriolis.cpp::relative_vorticity`。

## 4. 与 β 平面近似的关系

上一阶段的 `apply_coriolis` 采用 β 平面近似骨架（$f = 2\Omega\sin\phi$ 常数化），
忽略曲率项与动能梯度。本阶段升级为完整协变形式，含绝对涡度 $(f+\zeta)$
与动能梯度项，适用于全球模式。

## 5. 参考文献

- Putman, W. M., Lin, S.-J. (2007). *J. Comput. Phys.*, 227, 55–78.
- Sadourny, R. (1972). *Mon. Wea. Rev.*, 100, 136–144.
