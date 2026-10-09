// =============================================================================
//  科氏力与曲率项（Coriolis & Curvature，协变形式）
// =============================================================================
//  球面上协变形式的动量方程含度量 Christoffel 项（曲率项）。采用协变分量
//  存储时（Putman & Lin 2007 §3），动量方程的科氏力与曲率项可合并为：
//
//   ∂u_ξ/∂t  ←  (f + ζ) u_η - (∂Φ/∂ξ + ...)  ...   （协变形式）
//
//  其中 u_ξ、u_η 为协变风速分量，f = 2Ω sinφ 为科氏参数，ζ 为相对涡度。
//
//  关键几何项（gnomonic 投影，见 cubed_sphere_metric.md）：
//    - 相对涡度 ζ = (1/√G)(∂_ξ u_η - ∂_η u_ξ)（协变分量直接计算，无显式
//      Christoffel 项，Putman & Lin 2007 的核心优势）
//    - 动能梯度项 (1/2)∂_ξ(u·u)（用反变分量与度量）
//
//  本模块提供协变形式的科氏力 + 曲率项计算，替换 equations.cpp 中的
//  β 平面近似骨架。
//
//  参考文献：
//    Putman & Lin (2007), J. Comput. Phys. 227, 55–78（协变分量与曲率项）
//    Sadourny (1972), Mon. Wea. Rev. 100（球面守恒差分）
// =============================================================================

#pragma once

#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/grid/differential_operators.hpp"
#include "cubed_sph/common/types.hpp"

#include <vector>

namespace cubed_sph::dynamics {

// ---------------------------------------------------------------------------
// 科氏力与曲率项（协变形式）
// ---------------------------------------------------------------------------
class CoriolisCurvature {
public:
    CoriolisCurvature(const grid::CubedSphereGrid& grid, int nlev);

    // 计算科氏力 + 曲率项，累加到协变动量倾向 d(rho_u)/dt、d(rho_v)/dt。
    // 输入：
    //   u_cov, v_cov ：协变风速分量（由 ρu/ρ、ρv/ρ 还原）
    //   zeta          ：相对涡度（可预计算，或内部计算）
    // 输出累加：
    //   d_u, d_v      ：动量倾向（协变分量）
    void compute(const State& state,
                 const std::vector<StateReal>& u_cov,
                 const std::vector<StateReal>& v_cov,
                 const std::vector<StateReal>& zeta,
                 std::vector<StateReal>& d_u,
                 std::vector<StateReal>& d_v);

    // 相对涡度（协变分量直接计算，无显式 Christoffel 项）
    static void relative_vorticity(const grid::CubedSphereGrid& grid,
                                   int nlev,
                                   const std::vector<StateReal>& u_cov,
                                   const std::vector<StateReal>& v_cov,
                                   std::vector<StateReal>& zeta);

private:
    const grid::CubedSphereGrid& grid_;
    int nlev_;
};

}  // namespace cubed_sph::dynamics
