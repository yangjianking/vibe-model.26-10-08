// =============================================================================
//  半隐式 Helmholtz 算子（matrix-free 实现）
// =============================================================================
//  半隐式时间推进的核心：求解三维椭圆方程
//    (I - β² Δt² c_s² ∇²) δπ' = RHS
//
//  本模块实现 matrix-free 的离散 Helmholtz 算子 A(·)，以 stencil 回调形式
//  暴露给 BiCGStab 求解器（elliptic_solver.hpp）。算子结构：
//
//    A(δπ') = δπ' - β² Δt² [ ∇_h · (c_s² ∇_h δπ') + ∂_z(c_s² ∂_z δπ') ]
//
//  其中 ∇_h 为水平散度/梯度（立方球度量张量离散，见 differential_operators.hpp），
//  ∂_z 为垂直差分（地形跟随坐标的 ∂/∂z，经垂直 Jacobian 变换）。
//
//  离散形式（水平，立方球）：
//    ∇_h · (c_s² ∇_h φ) =
//      (1/√G) [ δ_ξ(√G c_s² g^{11} δ_ξ φ + √G c_s² g^{12} δ_η φ)
//              + δ_η(√G c_s² g^{21} δ_ξ φ + √G c_s² g^{22} δ_η φ) ]
//
//  垂直（地形跟随坐标）：
//    ∂_z(c_s² ∂_z φ) = (∂η/∂z) δ_η[ c_s² (∂η/∂z) δ_η φ ]
//
//  参考文献：
//    Wood et al. (2014), QJRMS 140, 1505–1528（ENDGame Helmholtz）
//    Buckeridge & Scheichl (2010), QJRMS 136, 2723–2735（球面多重网格）
// =============================================================================

#pragma once

#include "cubed_sph/grid/cubed_sphere_grid.hpp"
#include "cubed_sph/grid/vertical_coordinate.hpp"
#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/common/types.hpp"

#include <vector>

namespace cubed_sph::timeint {

// ---------------------------------------------------------------------------
// 半隐式 Helmholtz 算子（matrix-free）
// ---------------------------------------------------------------------------
// 缓存参考态声速平方 c_s²（与参考态位温 θ̄ 有关）与 β²Δt² 系数，提供
// A(x) = x - β²Δt² ∇·(c_s² ∇ x) 的应用。
class HelmholtzOperator {
public:
    HelmholtzOperator(const grid::CubedSphereGrid& grid,
                      const grid::VerticalCoordinate& vert,
                      int nlev,
                      StateReal beta_dt,          // β = α Δt
                      const std::vector<StateReal>& cs2);  // 参考态声速平方（每层）

    // 应用算子：out = A(x) = x - β²Δt² ∇·(c_s² ∇ x)
    // 向量布局与 State 一致：panel × n × n × nlev
    void apply(const std::vector<StateReal>& x, std::vector<StateReal>& out) const;

    // 提供对角近似（Jacobi 预条件用）：主对角 = 1 + β²Δt² * (水平+垂直系数)
    void diagonal(std::vector<StateReal>& diag) const;

    // 向量规模（与 State 字段规模一致）
    IIndex size() const;

private:
    const grid::CubedSphereGrid& grid_;
    const grid::VerticalCoordinate& vert_;
    int nlev_;
    StateReal beta_dt_;
    std::vector<StateReal> cs2_;   // 每层参考态声速平方

    // 水平散度算子（内部使用，作用于一层的场）
    void horizontal_laplacian(int k, const std::vector<StateReal>& phi,
                              std::vector<StateReal>& lap) const;
};

}  // namespace cubed_sph::timeint
