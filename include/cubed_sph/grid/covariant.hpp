// =============================================================================
//  协变/反变向量分量变换（Covariant/Contravariant Transform）
// =============================================================================
//  立方球 C 网格（Arakawa C/D 交错）上，风速分量以协变（covariant）或
//  反变（contravariant）分量存储。两种表示由度量张量相互转换：
//
//    u_contra^i = g^{ij} u_covar_j
//    u_covar_i  = g_{ij}  u_contra^j
//
//  其中 g_{ij} 为协变度量，g^{ij} 为反变度量（逆矩阵）。
//
//  采用协变分量存储的动机（Putman & Lin 2007）：
//   - 散度算子 ∇·u 在协变分量下形式简洁且严格守恒；
//   - 跨面板时协变分量按局部基向量旋转，旋转矩阵仅依赖面板朝向，
//     与位置无关（简化 halo 交换）。
//
//  参考文献：
//    Putman & Lin (2007), J. Comput. Phys. 227, 55–78（协变分量与旋转）
//    Sadourny (1972), Mon. Wea. Rev. 100（守恒散度离散）
// =============================================================================

#pragma once

#include "cubed_sph/grid/metric.hpp"

namespace cubed_sph::grid {

// ---------------------------------------------------------------------------
// 协变 → 反变： u_contra^i = g^{ij} u_covar_j
// ---------------------------------------------------------------------------
inline void covariant_to_contravariant(const MetricPoint& m,
                                       StateReal u_cov, StateReal v_cov,
                                       StateReal& u_contra, StateReal& v_contra) {
    u_contra = m.inv_g11 * u_cov + m.inv_g12 * v_cov;
    v_contra = m.inv_g12 * u_cov + m.inv_g22 * v_cov;
}

// ---------------------------------------------------------------------------
// 反变 → 协变： u_covar_i = g_{ij} u_contra^j
// ---------------------------------------------------------------------------
inline void contravariant_to_covariant(const MetricPoint& m,
                                       StateReal u_contra, StateReal v_contra,
                                       StateReal& u_cov, StateReal& v_cov) {
    u_cov = m.g11 * u_contra + m.g12 * v_contra;
    v_cov = m.g12 * u_contra + m.g22 * v_contra;
}

// ---------------------------------------------------------------------------
// 散度算子（守恒形式，协变分量输入）
// ---------------------------------------------------------------------------
// 在正交曲线坐标下，散度的守恒离散为：
//   ∇·u = (1/√G) [ ∂/∂ξ (√G u_contra^ξ) + ∂/∂η (√G u_contra^η) ]
// 调用方需先将协变分量转为反变分量，再按上式差分（见 dynamics/divergence）。
// 此处仅提供单点反变速度的辅助量 √G·u_contra^i（即"通量密度"），
// 供散度算子使用。
// ---------------------------------------------------------------------------
inline void flux_density(const MetricPoint& m,
                         StateReal u_contra, StateReal v_contra,
                         StateReal& sqrtG_u, StateReal& sqrtG_v) {
    sqrtG_u = m.sqrtG * u_contra;
    sqrtG_v = m.sqrtG * v_contra;
}

}  // namespace cubed_sph::grid
