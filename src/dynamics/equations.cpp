// =============================================================================
//  动力方程实现（显式右端项）
// =============================================================================

#include "cubed_sph/dynamics/equations.hpp"
#include "cubed_sph/dynamics/coriolis.hpp"
#include "cubed_sph/common/config.hpp"

#include <cmath>

namespace cubed_sph::dynamics {

using grid::kNumPanels;
using grid::MetricPoint;

Tendency::Tendency(const grid::CubedSphereGrid& g, int nlev)
    : drho(static_cast<size_t>(kNumPanels) * g.panel_n() * g.panel_n() * nlev, 0.0),
      drho_u(drho.size(), 0.0),
      drho_v(drho.size(), 0.0),
      drho_w(drho.size(), 0.0),
      drho_theta(drho.size(), 0.0) {}

void Tendency::zero() {
    std::fill(drho.begin(), drho.end(), 0.0);
    std::fill(drho_u.begin(), drho_u.end(), 0.0);
    std::fill(drho_v.begin(), drho_v.end(), 0.0);
    std::fill(drho_w.begin(), drho_w.end(), 0.0);
    std::fill(drho_theta.begin(), drho_theta.end(), 0.0);
}

DynamicsCore::DynamicsCore(const grid::CubedSphereGrid& grid,
                           const grid::VerticalCoordinate& vert,
                           int nlev,
                           const Config& cfg)
    : grid_(grid),
      vert_(vert),
      nlev_(nlev),
      coriolis_f_(cfg.get_double("dynamics.coriolis_f", 1.0e-4)) {}

// ---------------------------------------------------------------------------
// 显式右端项计算
// 说明：显式残差 = 水平气压梯度 + 科氏力/曲率 + 垂直方向（重力/垂直气压
//       梯度）。半隐式隐式项（快波线性项）由 timeint/ 在求解 Helmholtz 方程
//       后叠加。平流在 timeint 中单独调度（operator splitting）。
// ---------------------------------------------------------------------------
void DynamicsCore::compute_explicit(const State& state, Tendency& tend) {
    const int n = grid_.panel_n();
    const int nlev = nlev_;
    const StateReal dxi = 2.0 / static_cast<StateReal>(grid_.ncells());  // 等角坐标范围 [-1,1]

    tend.zero();

    // --- 水平气压梯度力（协变分量，Arakawa–Lamb 配对）---
    for (int p = 0; p < kNumPanels; ++p) {
        for (int k = 0; k < nlev; ++k) {
            for (int j = grid_.nhalo(); j < n - grid_.nhalo(); ++j) {
                for (int i = grid_.nhalo(); i < n - grid_.nhalo(); ++i) {
                    const IIndex idx = state.index(p, i, j, k);
                    const MetricPoint& m = grid_.metric(p, i, j);

                    const StateReal theta = state.rho_theta[idx] / state.rho[idx];
                    const StateReal dpi_dxi =
                        (state.exner[state.index(p, i + 1, j, k)] -
                         state.exner[state.index(p, i - 1, j, k)]) / (2.0 * dxi);
                    const StateReal dpi_deta =
                        (state.exner[state.index(p, i, j + 1, k)] -
                         state.exner[state.index(p, i, j - 1, k)]) / (2.0 * dxi);

                    // 协变气压梯度（乘逆度量缩放）
                    const StateReal pg_xi = -gas_constants::c_p * theta *
                                            (m.g11 * dpi_dxi + m.g12 * dpi_deta);
                    const StateReal pg_eta = -gas_constants::c_p * theta *
                                             (m.g12 * dpi_dxi + m.g22 * dpi_deta);

                    tend.drho_u[idx] += pg_xi;
                    tend.drho_v[idx] += pg_eta;
                }
            }
        }
    }

    // 科氏力与曲率项
    apply_coriolis(state, tend);

    // 垂直方向（重力、垂直气压梯度）
    apply_vertical(state, tend);
}

// ---------------------------------------------------------------------------
// 科氏力与曲率项（协变形式，完整实现见 coriolis.hpp/cpp）
// ---------------------------------------------------------------------------
// 采用协变分量计算，含度量 Christoffel 项（曲率项）与相对涡度
// （Putman & Lin 2007 §3）。本函数委托给 CoriolisCurvature 模块。
// ---------------------------------------------------------------------------
void DynamicsCore::apply_coriolis(const State& state, Tendency& tend) {
    const int nlev = nlev_;

    // 由动量密度还原协变风速 u_cov = ρu/ρ、v_cov = ρv/ρ
    std::vector<StateReal> u_cov(state.ncell(), 0.0);
    std::vector<StateReal> v_cov(state.ncell(), 0.0);
    for (IIndex i = 0; i < state.ncell(); ++i) {
        const StateReal rho = state.rho[i];
        u_cov[i] = (rho > 1e-12) ? state.rho_u[i] / rho : 0.0;
        v_cov[i] = (rho > 1e-12) ? state.rho_v[i] / rho : 0.0;
    }

    // 相对涡度（协变分量直接计算）
    std::vector<StateReal> zeta(state.ncell(), 0.0);
    CoriolisCurvature::relative_vorticity(grid_, nlev, u_cov, v_cov, zeta);

    // 科氏力 + 曲率项
    CoriolisCurvature cc(grid_, nlev);
    cc.compute(state, u_cov, v_cov, zeta, tend.drho_u, tend.drho_v);
}

// ---------------------------------------------------------------------------
// 垂直方向项（重力 + 垂直气压梯度）
// ---------------------------------------------------------------------------
// 垂直动量方程（非静力，守恒形式）：
//   ∂(ρw)/∂t = -∂p/∂z - ρg + ...（垂直平流在平流模块单独处理）
//
// 用 Exner 压力 π 与位温 θ 表示垂直气压梯度（Harris et al. 2021 §4）：
//   -∂p/∂z = -c_p ρ θ ∂π/∂z - ρ g
//
// 采用扰动形式（在参考态附近）：
//   -∂p'/∂z - ρ' g ≈ -c_p (ρθ)' ∂π̄/∂z - c_p ρ̄ θ̄ ∂π'/∂z - ρ' g
//
// 本函数计算显式的垂直气压梯度项（用 Exner 压力与位温的离散形式）与重力项。
// 半隐式隐式部分（∂π'/∂z 与浮力项的耦合）由 Helmholtz 求解器处理。
// 参考文献：Wood et al. (2014, QJRMS 140)；Harris et al. (2021, GFDL TM §5)。
// ---------------------------------------------------------------------------
void DynamicsCore::apply_vertical(const State& state, Tendency& tend) {
    const int n = grid_.panel_n();
    const int nlev = nlev_;

    for (int p = 0; p < kNumPanels; ++p) {
        for (int j = grid_.nhalo(); j < n - grid_.nhalo(); ++j) {
            for (int i = grid_.nhalo(); i < n - grid_.nhalo(); ++i) {
                for (int k = 0; k < nlev; ++k) {
                    const IIndex idx = state.index(p, i, j, k);
                    const StateReal rho = state.rho[idx];
                    const StateReal theta = state.rho_theta[idx] / rho;

                    // 垂直 Exner 压力梯度 ∂π/∂z
                    // 用相邻层中心差分（边界层单侧差分）
                    StateReal dpi_dz = 0.0;
                    if (k == 0) {
                        const StateReal dz = vert_.z_full(1) - vert_.z_full(0);
                        dpi_dz = (state.exner[state.index(p, i, j, 1)] -
                                  state.exner[idx]) / dz;
                    } else if (k == nlev - 1) {
                        const StateReal dz = vert_.z_full(nlev - 1) -
                                             vert_.z_full(nlev - 2);
                        dpi_dz = (state.exner[idx] -
                                  state.exner[state.index(p, i, j, nlev - 2)]) / dz;
                    } else {
                        const StateReal dz = vert_.z_full(k + 1) -
                                             vert_.z_full(k - 1);
                        dpi_dz = (state.exner[state.index(p, i, j, k + 1)] -
                                  state.exner[state.index(p, i, j, k - 1)]) / dz;
                    }

                    // 垂直气压梯度力（-c_p ρ θ ∂π/∂z）与重力（-ρg）
                    const StateReal pg_z = -gas_constants::c_p * rho * theta * dpi_dz;
                    const StateReal gravity = -rho * gas_constants::g;

                    tend.drho_w[idx] += pg_z + gravity;
                }
            }
        }
    }
}

}  // namespace cubed_sph::dynamics
