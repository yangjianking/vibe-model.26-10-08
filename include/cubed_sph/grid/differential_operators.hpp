// =============================================================================
//  立方球网格差分算子（梯度 / 散度 / 涡度 / 拉普拉斯）
// =============================================================================
//  在 Arakawa C 交错网格上定义水平差分算子。所有算子均以守恒/配对的
//  离散形式给出，保证散度与压力梯度的能量配对（Arakawa–Lamb 1981）。
//
//  离散约定（单元中心标量 φ，面通量 u_contra）：
//   散度（守恒）：  ∇·u = (1/√G) [ δ_ξ(√G·u^ξ) + δ_η(√G·u^η) ]
//   梯度（标量）：  (∂φ/∂ξ, ∂φ/∂η) 在面/中心交错处用中心差分
//   拉普拉斯：      ∇²φ = ∇·(∇φ)，二阶（可选四阶）
//
//  参考文献：
//    Sadourny (1972), Mon. Wea. Rev. 100
//    Arakawa & Lamb (1981), Mon. Wea. Rev. 109, 18–36（能量守恒配对）
//    Putman & Lin (2007), J. Comput. Phys. 227, 55–78
// =============================================================================

#pragma once

#include "cubed_sph/grid/cubed_sphere_grid.hpp"

#include <vector>

namespace cubed_sph::grid {

// ---------------------------------------------------------------------------
// 差分算子结果容器（与网格同布局：panel × n × n）
// ---------------------------------------------------------------------------
using Field2D = std::vector<StateReal>;

// ---------------------------------------------------------------------------
// 二阶中心差分模板（一维，周期/面板内有效，面板边界由 halo 保证连续性）
// ---------------------------------------------------------------------------
// δ_ξ f 在 i+1/2 处（面）： (f(i+1) - f(i)) / Δξ
inline StateReal ddx_face(const CubedSphereGrid& g, const Field2D& f,
                          int panel, int i, int j) {
    const int n = g.panel_n();
    const StateReal dxi = 2.0 / static_cast<StateReal>(g.ncells());  // 等角坐标范围 [-1,1]
    const auto idx = [&](int p, int ii, int jj) {
        return static_cast<size_t>(p) * n * n + jj * n + ii;
    };
    return (f[idx(panel, i + 1, j)] - f[idx(panel, i, j)]) / dxi;
}

// ---------------------------------------------------------------------------
// 散度（守恒形式）
// ---------------------------------------------------------------------------
// 输入：u_contra, v_contra 为反变速度分量（已乘 √G 与否由 use_flux 控制）。
// 本函数假设输入为"通量密度" √G·u^ξ、√G·v^η（由 covariant.hpp 的
// flux_density 得到），输出散度到单元中心。
//
//   (∇·u)_{i,j} = (1/√G_{i,j}) [ (√G u^ξ)_{i+1/2} - (√G u^ξ)_{i-1/2}
//                               + (√G v^η)_{j+1/2} - (√G v^η)_{j-1/2} ] / Δξ
inline void divergence(const CubedSphereGrid& g,
                       const std::vector<StateReal>& sqrtG_u,   // 面通量 ξ
                       const std::vector<StateReal>& sqrtG_v,   // 面通量 η
                       int panel,
                       Field2D& div) {
    const int n = g.panel_n();
    const StateReal dxi = 2.0 / static_cast<StateReal>(g.ncells());  // 等角坐标范围 [-1,1]
    const auto idx = [&](int p, int ii, int jj) {
        return static_cast<size_t>(p) * n * n + jj * n + ii;
    };
    for (int j = g.nhalo(); j < n - g.nhalo(); ++j) {
        for (int i = g.nhalo(); i < n - g.nhalo(); ++i) {
            const MetricPoint& m = g.metric(panel, i, j);
            const StateReal du = sqrtG_u[idx(panel, i + 1, j)] -
                                 sqrtG_u[idx(panel, i, j)];
            const StateReal dv = sqrtG_v[idx(panel, i, j + 1)] -
                                 sqrtG_v[idx(panel, i, j)];
            div[idx(panel, i, j)] = (du + dv) / (m.sqrtG * dxi);
        }
    }
}

// ---------------------------------------------------------------------------
// 标量梯度（中心 → 面）
// ---------------------------------------------------------------------------
// 在 C 网格上，压力梯度力作用于面（与速度共位），由中心压力差分得到。
inline void gradient_center_to_face(const CubedSphereGrid& g,
                                    const Field2D& phi, int panel,
                                    std::vector<StateReal>& grad_xi,  // 东面
                                    std::vector<StateReal>& grad_eta) { // 北面
    const int n = g.panel_n();
    const StateReal dxi = 2.0 / static_cast<StateReal>(g.ncells());  // 等角坐标范围 [-1,1]
    const auto idx = [&](int p, int ii, int jj) {
        return static_cast<size_t>(p) * n * n + jj * n + ii;
    };
    for (int j = g.nhalo(); j < n - g.nhalo(); ++j) {
        for (int i = g.nhalo(); i < n - g.nhalo(); ++i) {
            grad_xi[idx(panel, i, j)] = (phi[idx(panel, i, j)] -
                                         phi[idx(panel, i - 1, j)]) / dxi;
            grad_eta[idx(panel, i, j)] = (phi[idx(panel, i, j)] -
                                          phi[idx(panel, i, j - 1)]) / dxi;
        }
    }
}

// ---------------------------------------------------------------------------
// 相对涡度（C 网格 → 顶点/D 网格）
// ---------------------------------------------------------------------------
// 涡度 ζ = (1/√G)(δ_ξ v_cov - δ_η u_cov)，用协变分量直接计算可避免度量
// 的显式 Christoffel 项（Putman & Lin 2007）。
inline void vorticity(const CubedSphereGrid& g,
                      const std::vector<StateReal>& u_cov,
                      const std::vector<StateReal>& v_cov,
                      int panel, Field2D& zeta) {
    const int n = g.panel_n();
    const StateReal dxi = 2.0 / static_cast<StateReal>(g.ncells());  // 等角坐标范围 [-1,1]
    const auto idx = [&](int p, int ii, int jj) {
        return static_cast<size_t>(p) * n * n + jj * n + ii;
    };
    for (int j = g.nhalo(); j < n - g.nhalo(); ++j) {
        for (int i = g.nhalo(); i < n - g.nhalo(); ++i) {
            const MetricPoint& m = g.metric(panel, i, j);
            const StateReal dv_dxi = (v_cov[idx(panel, i + 1, j)] -
                                      v_cov[idx(panel, i, j)]) / dxi;
            const StateReal du_deta = (u_cov[idx(panel, i, j + 1)] -
                                       u_cov[idx(panel, i, j)]) / dxi;
            zeta[idx(panel, i, j)] = (dv_dxi - du_deta) / m.sqrtG;
        }
    }
}

}  // namespace cubed_sph::grid
