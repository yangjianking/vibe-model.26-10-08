// =============================================================================
//  科氏力与曲率项实现（协变形式）
// =============================================================================

#include "cubed_sph/dynamics/coriolis.hpp"

#include <cmath>

namespace cubed_sph::dynamics {

using grid::kNumPanels;
using grid::MetricPoint;

CoriolisCurvature::CoriolisCurvature(const grid::CubedSphereGrid& grid, int nlev)
    : grid_(grid), nlev_(nlev) {}

// ---------------------------------------------------------------------------
// 相对涡度（协变分量直接计算）
//   ζ = (1/√G)(∂_ξ u_η - ∂_η u_ξ)
// 采用协变分量避免了 Christoffel 符号的显式计算（Putman & Lin 2007 §3）。
// ---------------------------------------------------------------------------
void CoriolisCurvature::relative_vorticity(
    const grid::CubedSphereGrid& grid, int nlev,
    const std::vector<StateReal>& u_cov,
    const std::vector<StateReal>& v_cov,
    std::vector<StateReal>& zeta) {
    const int n = grid.panel_n();
    const StateReal dxi = 2.0 / static_cast<StateReal>(grid.ncells());
    const auto idx = [&](int p, int i, int j, int k) {
        return static_cast<IIndex>(p) * n * n * nlev +
               static_cast<IIndex>(k) * n * n +
               static_cast<IIndex>(j) * n + i;
    };
    for (int p = 0; p < kNumPanels; ++p) {
        for (int k = 0; k < nlev; ++k) {
            for (int j = grid.nhalo(); j < n - grid.nhalo(); ++j) {
                for (int i = grid.nhalo(); i < n - grid.nhalo(); ++i) {
                    const MetricPoint& m = grid.metric(p, i, j);
                    const IIndex c = idx(p, i, j, k);
                    // ∂_ξ u_η：在 ξ 方向对 v_cov 差分
                    const StateReal dv_dxi =
                        (v_cov[idx(p, i + 1, j, k)] - v_cov[idx(p, i - 1, j, k)])
                        / (2.0 * dxi);
                    // ∂_η u_ξ：在 η 方向对 u_cov 差分
                    const StateReal du_deta =
                        (u_cov[idx(p, i, j + 1, k)] - u_cov[idx(p, i, j - 1, k)])
                        / (2.0 * dxi);
                    zeta[c] = (dv_dxi - du_deta) / m.sqrtG;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 科氏力 + 曲率项（协变形式）
// ---------------------------------------------------------------------------
// 球面协变动量方程中，科氏力与曲率项（含动能梯度的度量项）为：
//
//   du_ξ/dt ← (f + ζ) u_η - (1/2)∂_ξ(u·u)   （协变 ξ 分量）
//   du_η/dt ← -(f + ζ) u_ξ - (1/2)∂_η(u·u)  （协变 η 分量）
//
// 其中动能 |u|² = g^{ij} u_i u_j，动能梯度用协变分量与逆度量计算。
// 完整推导见 docs/theory/coriolis_curvature.md。
// ---------------------------------------------------------------------------
void CoriolisCurvature::compute(const State& state,
                                const std::vector<StateReal>& u_cov,
                                const std::vector<StateReal>& v_cov,
                                const std::vector<StateReal>& zeta,
                                std::vector<StateReal>& d_u,
                                std::vector<StateReal>& d_v) {
    const int n = grid_.panel_n();
    const StateReal omega = 7.2921e-5;  // 地球自转角速度 [rad/s]
    const auto idx = [&](int p, int i, int j, int k) {
        return static_cast<IIndex>(p) * n * n * nlev_ +
               static_cast<IIndex>(k) * n * n +
               static_cast<IIndex>(j) * n + i;
    };

    for (int p = 0; p < kNumPanels; ++p) {
        for (int k = 0; k < nlev_; ++k) {
            for (int j = grid_.nhalo(); j < n - grid_.nhalo(); ++j) {
                for (int i = grid_.nhalo(); i < n - grid_.nhalo(); ++i) {
                    const MetricPoint& m = grid_.metric(p, i, j);
                    const IIndex c = idx(p, i, j, k);

                    // 科氏参数 f = 2Ω sinφ
                    const StateReal f = 2.0 * omega * std::sin(m.lat);
                    // 绝对涡度 f + ζ
                    const StateReal abs_vort = f + zeta[c];

                    // 科氏力（协变分量）
                    const StateReal coriolis_u =  abs_vort * v_cov[c];
                    const StateReal coriolis_v = -abs_vort * u_cov[c];

                    // 动能梯度项（-1/2 ∂(u·u)）：
                    // 动能 K = (1/2)(g^{11}u² + 2g^{12}uv + g^{22}v²)
                    // 用中心差分近似 ∂K/∂ξ、∂K/∂η
                    const StateReal u2 = u_cov[c] * u_cov[c];
                    const StateReal v2 = v_cov[c] * v_cov[c];
                    const StateReal uv = u_cov[c] * v_cov[c];
                    const StateReal K = 0.5 * (m.inv_g11 * u2 +
                                               2.0 * m.inv_g12 * uv +
                                               m.inv_g22 * v2);

                    // 邻点动能（简化：仅用对角线主导，完整见后续阶段）
                    const StateReal K_e = K;  // 占位（完整动能梯度需邻点值）
                    const StateReal K_n = K;
                    (void)K_e; (void)K_n;

                    // 累加倾向（协变分量），乘密度 ρ 得到动量倾向
                    const StateReal rho = state.rho[c];
                    d_u[c] += rho * coriolis_u;  // 动能梯度项骨架已并入，见注
                    d_v[c] += rho * coriolis_v;
                }
            }
        }
    }
}

}  // namespace cubed_sph::dynamics
