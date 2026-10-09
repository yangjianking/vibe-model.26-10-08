// =============================================================================
//  立方球网格度量（Metric）：gnomonic equiangular 投影的几何量
// =============================================================================
//  立方球的核心几何：球面 → 立方体面的 gnomonic（等角球心）投影。
//
//  设面板局部坐标 (ξ, η) ∈ [-1, 1]²（等角坐标，ξ = tan θ，θ ∈ [-π/4, π/4]），
//  面板中心对应 ξ = η = 0。球面上一点的笛卡尔坐标由投影反算得到。
//
//  关键几何量（本结构体缓存，供差分算子复用）：
//    - g_ij   ：度量张量（协变分量），用于协变/反变变换与散度/梯度
//    - sqrtG  ：Jacobian 行列式平方根 √G = sqrt(det(g_ij))，即面积元缩放
//    - inv_g_ij：逆度量张量（反变分量）
//    - christoffel：Christoffel 符号（可选，用于协变导数与向量传输）
//
//  推导细节见 docs/theory/cubed_sphere_metric.md，此处仅缓存数值结果。
//
//  参考文献：
//    Sadourny (1972), Mon. Wea. Rev. 100, 136–144（守恒有限差分，球面）
//    Ronchi et al. (1996), J. Comput. Phys. 124（gnomonic 投影度量）
//    Putman & Lin (2007), J. Comput. Phys. 227, 55–78（度量与交错）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"

namespace cubed_sph::grid {

// ---------------------------------------------------------------------------
// 单个网格点（或单元中心/面/边）的几何度量
// ---------------------------------------------------------------------------
struct MetricPoint {
    // 等角坐标 (ξ, η)，无量纲，范围 [-1, 1]
    // 定义：ξ = tan(θ_ξ)，其中 θ_ξ ∈ [-π/4, π/4] 为面板中心到该点在 e1 方向
    // 的球面角。故 ξ ∈ [-1, 1]（gnomonic equiangular 投影，Ronchi et al. 1996）。
    StateReal xi  = 0.0;
    StateReal eta = 0.0;

    // 度量张量（协变分量），2x2 对称
    StateReal g11 = 1.0;
    StateReal g12 = 0.0;
    StateReal g22 = 1.0;

    // 逆度量张量（反变分量），2x2 对称
    StateReal inv_g11 = 1.0;
    StateReal inv_g12 = 0.0;
    StateReal inv_g22 = 1.0;

    // Jacobian 行列式平方根 √G（面积元 dA = √G dξ dη）
    StateReal sqrtG = 1.0;

    // 球面单位法向（面板局部，指向球外），用于协变↔反变与向量旋转
    StateReal normal_x = 0.0;
    StateReal normal_y = 0.0;
    StateReal normal_z = 1.0;

    // 对应球面点的笛卡尔单位坐标（半径归一化到 1）
    StateReal cart_x = 1.0;
    StateReal cart_y = 0.0;
    StateReal cart_z = 0.0;

    // 经纬度（弧度），便于诊断与 I/O 的 CF 约定输出
    StateReal lon = 0.0;
    StateReal lat = 0.0;
};

// ---------------------------------------------------------------------------
// 计算单个 (ξ, η) 点的完整度量（gnomonic equiangular 投影）
// ---------------------------------------------------------------------------
// 投影公式（Ronchi et al. 1996, Eq. 4-8）：
//   令 r = sqrt(1 + ξ² + η²)，面板中心方向为面板法向 n̂。
//   球面单位向量： x = (1, ξ, η)/r  （面板局部正交基下的分量）
//   经纬度由面板朝向矩阵旋转到地理系（见 panel_orientation）。
//
// 注意：gnomonic equiangular 投影的等角坐标 ξ = tan(θ_ξ) 为无量纲量，
//   范围 [-1, 1]（θ_ξ ∈ [-π/4, π/4]）。度量张量（Sadourny 1972；
//   Putman & Lin 2007, Eq. A1-A6）：
//     g11 = (1 + η²) / r⁴,  g22 = (1 + ξ²) / r⁴,  g12 = -ξ η / r⁴
//   其中 r² = 1 + ξ² + η²。
//   √G = 1 / r³
//
// 详细推导见 docs/theory/cubed_sphere_metric.md。
// ---------------------------------------------------------------------------
MetricPoint compute_metric(StateReal xi, StateReal eta,
                           int panel_orientation_index);

}  // namespace cubed_sph::grid
