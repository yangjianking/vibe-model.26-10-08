# 垂直坐标与地形跟随（Vertical Coordinate）

> 模块：`src/grid/` · 对应代码：`vertical_coordinate.hpp/cpp`、`topography.hpp/cpp`

## 1. 问题陈述

垂直方向需在模式顶（约 30 km）与地形表面之间离散。纯高度坐标在地形处边界
复杂（shaved-cell）；地形跟随坐标将地形表面映射为坐标面，边界条件简单，但
要求坐标面单调不交叉。

## 2. 地形跟随坐标（Schär et al. 2002 平滑形式）

归一化垂直坐标 $\eta \in [0, 1]$（0 = 地面，1 = 模式顶），高度：

$$
z(\eta) = z_{\text{top}} \eta + h_s \, b(\eta)
$$

其中 $h_s$ 为平滑地形高度，$b(\eta)$ 为地形衰减函数：

$$
b(\eta) = \frac{\sinh\!\left(\frac{\eta_{\text{top}} - \eta}{s}\right)}
             {\sinh\!\left(\frac{\eta_{\text{top}}}{s}\right)}
$$

满足 $b(0) = 1$（地面）、$b(\eta_{\text{top}}) = 0$（模式顶）。参数 $s$ 控制
地形衰减尺度：$s$ 越小，地形影响越局限于近地面，高层坐标面越接近纯高度坐标。

## 3. 度量项（垂直 Jacobian）

$$
J = \frac{\partial z}{\partial \eta}
= z_{\text{top}} + h_s \frac{db}{d\eta}
$$

其中

$$
\frac{db}{d\eta} = -\frac{1}{s}\frac{\cosh\!\left(\frac{\eta_{\text{top}}-\eta}{s}\right)}
                             {\sinh\!\left(\frac{\eta_{\text{top}}}{s}\right)}
$$

垂直导数的坐标变换：$\partial/\partial z = (1/J)\,\partial/\partial\eta$。

## 4. 坐标面单调性条件

坐标面不交叉要求 $J > 0$，即

$$
z_{\text{top}} + h_s \frac{db}{d\eta} > 0
\qquad \Longleftrightarrow \qquad
h_s < \frac{z_{\text{top}}\, s\, \sinh(1/s)}{\cosh(1/s)}
$$

对典型参数 $z_{\text{top}} = 30000$ m、$s = 5$，该上界约为 $30000 \times 5
\times \tanh(0.2) \approx 29600$ m，远超地球最大地形（8848 m），故单调性
恒成立（已数值验证，见 `vertical_coordinate.cpp` 配套验证）。

## 5. 地形斜率限制（Janjić 1984 类判据）

限制地形坡度 $|\partial h/\partial x| \le s_{\max}$（可配，默认 0.5），超出
部分削峰，避免陡峭地形引起虚假的地形强迫与数值不稳定。实现见
`topography.hpp::limit_slope`。

## 6. 交错方案

- **Lorenz**：速度、温度、气压同位（经典，实现简单）。
- **Charney–Phillips**：温度与垂直速度交错（避免虚假静力不平衡，本框架默认）。

## 7. 参考文献

- Phillips, N. A. (1957). A coordinate system having some special advantages
  for numerical forecasting. *J. Meteorol.*, 14, 184–185.
- Simmons, A. J., Burridge, D. M. (1981). An energy and angular-momentum
  conserving vertical finite-difference scheme and hybrid vertical coordinates.
  *Mon. Wea. Rev.*, 109, 758–766.
- Schär, C., et al. (2002). A new terrain-following vertical coordinate
  formulation for atmospheric prediction models. *Mon. Wea. Rev.*, 130, 2459–2480.
- Janjić, Z. I. (1984). [TO-VERIFY]（地形斜率限制，待补充确切出处）。
