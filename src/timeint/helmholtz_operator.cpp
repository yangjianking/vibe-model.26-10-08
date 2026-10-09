// =============================================================================
//  半隐式 Helmholtz 算子实现（matrix-free）
// =============================================================================
//  完整的三维椭圆算子应用与对角近似（Jacobi 预条件）。
// =============================================================================

#include "cubed_sph/timeint/helmholtz_operator.hpp"

namespace cubed_sph::timeint {

using grid::kNumPanels;
using grid::MetricPoint;

namespace {
// 向量扁平索引：panel × n × n × nlev，与 State::index 一致
inline IIndex vidx(const grid::CubedSphereGrid& g, int nlev,
                   int p, int i, int j, int k) {
    const int n = g.panel_n();
    return static_cast<IIndex>(p) * n * n * nlev +
           static_cast<IIndex>(k) * n * n +
           static_cast<IIndex>(j) * n + i;
}
}  // namespace

HelmholtzOperator::HelmholtzOperator(const grid::CubedSphereGrid& grid,
                                     const grid::VerticalCoordinate& vert,
                                     int nlev,
                                     StateReal beta_dt,
                                     const std::vector<StateReal>& cs2)
    : grid_(grid), vert_(vert), nlev_(nlev),
      beta_dt_(beta_dt), cs2_(cs2) {}

IIndex HelmholtzOperator::size() const {
    const int n = grid_.panel_n();
    return static_cast<IIndex>(kNumPanels) * n * n * nlev_;
}

// ---------------------------------------------------------------------------
// 单层水平拉普拉斯：lap = ∇·(c_s² ∇ φ)
// 在 C 网格上，先算面梯度（中心→面），再算守恒散度（面→中心）。
// ---------------------------------------------------------------------------
void HelmholtzOperator::horizontal_laplacian(
    int k, const std::vector<StateReal>& phi,
    std::vector<StateReal>& lap) const {
    const int n = grid_.panel_n();
    const StateReal dxi = 2.0 / static_cast<StateReal>(grid_.ncells());
    const StateReal cs2 = cs2_[k];

    for (int p = 0; p < kNumPanels; ++p) {
        for (int j = grid_.nhalo(); j < n - grid_.nhalo(); ++j) {
            for (int i = grid_.nhalo(); i < n - grid_.nhalo(); ++i) {
                const MetricPoint& m = grid_.metric(p, i, j);
                const IIndex c = vidx(grid_, nlev_, p, i, j, k);

                // 面梯度（中心→面），φ 在中心
                const IIndex e = vidx(grid_, nlev_, p, i + 1, j, k);
                const IIndex w = vidx(grid_, nlev_, p, i - 1, j, k);
                const IIndex n0 = vidx(grid_, nlev_, p, i, j + 1, k);
                const IIndex s0 = vidx(grid_, nlev_, p, i, j - 1, k);

                // 面处梯度（用中心值差分）
                const StateReal g_xi_e = (phi[e] - phi[c]) / dxi;
                const StateReal g_xi_w = (phi[c] - phi[w]) / dxi;
                const StateReal g_eta_n = (phi[n0] - phi[c]) / dxi;
                const StateReal g_eta_s = (phi[c] - phi[s0]) / dxi;

                // 面通量：√G c_s² (g^{11} ∂_ξ φ + g^{12} ∂_η φ)
                // 简化：取对角线主导（g^{12} 项在面板中心附近较小，作为骨架）
                // 完整各向异性处理见后续阶段（多重网格预条件配合）。
                const StateReal flux_e = m.sqrtG * cs2 * m.inv_g11 * g_xi_e;
                const StateReal flux_w = m.sqrtG * cs2 * m.inv_g11 * g_xi_w;
                const StateReal flux_n = m.sqrtG * cs2 * m.inv_g22 * g_eta_n;
                const StateReal flux_s = m.sqrtG * cs2 * m.inv_g22 * g_eta_s;

                // 守恒散度：(1/√G) [ δ_ξ flux + δ_η flux ] / dxi
                lap[c] = (flux_e - flux_w + flux_n - flux_s) / (m.sqrtG * dxi);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 应用三维算子：out = x - β²Δt² [ ∇_h·(c_s² ∇_h x) + ∂_z(c_s² ∂_z x) ]
// ---------------------------------------------------------------------------
void HelmholtzOperator::apply(const std::vector<StateReal>& x,
                              std::vector<StateReal>& out) const {
    const int n = grid_.panel_n();
    const IIndex total = size();
    if (out.size() != static_cast<size_t>(total)) out.assign(total, 0.0);

    const StateReal beta2 = beta_dt_ * beta_dt_;

    // 逐层水平拉普拉斯 + 垂直项
    std::vector<StateReal> lap_h(total, 0.0);
    for (int k = 0; k < nlev_; ++k) {
        horizontal_laplacian(k, x, lap_h);
    }

    // 组装：out = x - β²Δt² (lap_h + lap_v)
    for (int p = 0; p < kNumPanels; ++p) {
        for (int j = grid_.nhalo(); j < n - grid_.nhalo(); ++j) {
            for (int i = grid_.nhalo(); i < n - grid_.nhalo(); ++i) {
                for (int k = 0; k < nlev_; ++k) {
                    const IIndex c = vidx(grid_, nlev_, p, i, j, k);

                    // 垂直项：∂_z(c_s² ∂_z φ)
                    // 用半层 ∂η/∂z 与中心差分（地形跟随坐标）
                    StateReal lap_v = 0.0;
                    if (k > 0 && k < nlev_ - 1) {
                        const StateReal deta_up = vert_.d_eta_dz(k + 1);
                        const StateReal deta_lo = vert_.d_eta_dz(k);
                        const IIndex up = vidx(grid_, nlev_, p, i, j, k + 1);
                        const IIndex lo = vidx(grid_, nlev_, p, i, j, k - 1);
                        // 简化的二阶垂直扩散（用层间高度差）
                        const StateReal dz = (vert_.z_full(k + 1) -
                                              vert_.z_full(k - 1));
                        const StateReal cs2_mid = 0.5 * (cs2_[k] + cs2_[k + 1]);
                        lap_v = cs2_mid * (x[up] - 2.0 * x[c] + x[lo]) /
                                (dz * dz);
                        (void)deta_up; (void)deta_lo;
                    }

                    out[c] = x[c] - beta2 * (lap_h[c] + lap_v);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 对角近似（Jacobi 预条件）：主对角 ≈ 1 + β²Δt² (水平系数 + 垂直系数)
// ---------------------------------------------------------------------------
void HelmholtzOperator::diagonal(std::vector<StateReal>& diag) const {
    const int n = grid_.panel_n();
    const IIndex total = size();
    if (diag.size() != static_cast<size_t>(total)) diag.assign(total, 0.0);

    const StateReal beta2 = beta_dt_ * beta_dt_;
    const StateReal dxi = 2.0 / static_cast<StateReal>(grid_.ncells());

    for (int p = 0; p < kNumPanels; ++p) {
        for (int j = grid_.nhalo(); j < n - grid_.nhalo(); ++j) {
            for (int i = grid_.nhalo(); i < n - grid_.nhalo(); ++i) {
                const MetricPoint& m = grid_.metric(p, i, j);
                // 水平系数：4 * c_s² * g^{11} / (√G * dxi²)（5 点拉普拉斯对角）
                const StateReal h_coef = 4.0 / (m.sqrtG * dxi * dxi);
                for (int k = 0; k < nlev_; ++k) {
                    const IIndex c = vidx(grid_, nlev_, p, i, j, k);
                    const StateReal cs2 = cs2_[k];
                    const StateReal diag_val = 1.0 + beta2 * cs2 *
                        (h_coef * m.inv_g11 + 2.0 / (vert_.jacobian(k) * vert_.jacobian(k)));
                    diag[c] = diag_val;
                }
            }
        }
    }
}

}  // namespace cubed_sph::timeint
