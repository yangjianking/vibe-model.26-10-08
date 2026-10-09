// =============================================================================
//  守恒诊断实现
// =============================================================================

#include "cubed_sph/diagnostics/diagnostics.hpp"
#include "cubed_sph/dynamics/thermo.hpp"

#include <cmath>
#include <algorithm>

namespace cubed_sph::diagnostics {

using grid::kNumPanels;
using grid::MetricPoint;
namespace gas_constants = dynamics::gas_constants;

// ---------------------------------------------------------------------------
// 全球积分（长双精度补偿求和 + gnomonic 度量加权）
// ---------------------------------------------------------------------------
StateReal global_integral(const grid::CubedSphereGrid& grid,
                          const grid::VerticalCoordinate& vert,
                          const std::vector<StateReal>& field) {
    const int n = grid.panel_n();
    const int nlev = vert.nlev();
    const int nhalo = grid.nhalo();
    const StateReal dxi = 2.0 / static_cast<StateReal>(grid.ncells());
    const StateReal R2 = grid.radius() * grid.radius();  // 单位球 → 物理球面积元因子

    long double acc = 0.0L;
    for (int p = 0; p < kNumPanels; ++p) {
        for (int k = 0; k < nlev; ++k) {
            const StateReal dz = vert.z_half(k + 1) - vert.z_half(k);  // 层厚
            for (int j = nhalo; j < n - nhalo; ++j) {
                for (int i = nhalo; i < n - nhalo; ++i) {
                    const MetricPoint& m = grid.metric(p, i, j);
                    const size_t idx = static_cast<size_t>(p) * n * n * nlev +
                                       static_cast<size_t>(k) * n * n +
                                       static_cast<size_t>(j) * n + i;
                    // 积分元：R² · √G dξ dη dz（等角坐标 dξ=dη=dxi）
                    acc += (long double)field[idx] * (long double)m.sqrtG *
                           (long double)(R2 * dxi * dxi) * (long double)dz;
                }
            }
        }
    }
    return static_cast<StateReal>(acc);
}

// ---------------------------------------------------------------------------
// ConservationDiagnostics
// ---------------------------------------------------------------------------
ConservationDiagnostics::ConservationDiagnostics(
    const grid::CubedSphereGrid& grid,
    const grid::VerticalCoordinate& vert)
    : grid_(grid), vert_(vert) {}

ConservationSnapshot ConservationDiagnostics::compute(
    const dynamics::State& state) const {
    ConservationSnapshot snap;

    const int n = grid_.panel_n();
    const int nlev = vert_.nlev();
    const int nhalo = grid_.nhalo();
    const StateReal dxi = 2.0 / static_cast<StateReal>(grid_.ncells());
    const StateReal radius = grid_.radius();
    const StateReal R2 = radius * radius;  // 单位球 → 物理球面积元因子
    const StateReal omega = 7.2921150e-5;  // 地球自转角速度 [rad/s]

    long double mass = 0.0L;
    long double theta_int = 0.0L;
    long double ke = 0.0L;
    long double ie = 0.0L;
    long double pe = 0.0L;
    long double angmom = 0.0L;

    for (int p = 0; p < kNumPanels; ++p) {
        for (int k = 0; k < nlev; ++k) {
            const StateReal dz = vert_.z_half(k + 1) - vert_.z_half(k);
            const StateReal z = vert_.z_full(k);
            for (int j = nhalo; j < n - nhalo; ++j) {
                for (int i = nhalo; i < n - nhalo; ++i) {
                    const MetricPoint& m = grid_.metric(p, i, j);
                    const size_t idx = static_cast<size_t>(p) * n * n * nlev +
                                       static_cast<size_t>(k) * n * n +
                                       static_cast<size_t>(j) * n + i;

                    const StateReal rho = state.rho[idx];
                    const StateReal rho_theta = state.rho_theta[idx];
                    const StateReal rho_u = state.rho_u[idx];
                    const StateReal rho_v = state.rho_v[idx];
                    const StateReal rho_w = state.rho_w[idx];

                    const StateReal theta = (rho > 1e-12)
                        ? rho_theta / rho : 0.0;
                    // 温度 T = θ π（π 由 Exner 或诊断；此处用 state.exner）
                    const StateReal T = theta * state.exner[idx];
                    // 位势 Φ = g z（state.geopotential 或直接 gz）
                    const StateReal Phi = gas_constants::g * z;

                    // 速度（协变 → 反变需度量；动能用协变分量近似：
                    //  |v|² ≈ u_cov·u_contra = u_i g^{ij} u_j）
                    // 这里用逆度量还原反变速度，再算 |v|²。
                    const StateReal u_cov = (rho > 1e-12) ? rho_u / rho : 0.0;
                    const StateReal v_cov = (rho > 1e-12) ? rho_v / rho : 0.0;
                    const StateReal w_cov = (rho > 1e-12) ? rho_w / rho : 0.0;
                    const StateReal u_contra = m.inv_g11 * u_cov + m.inv_g12 * v_cov;
                    const StateReal v_contra = m.inv_g12 * u_cov + m.inv_g22 * v_cov;
                    const StateReal v2 = u_cov * u_contra + v_cov * v_contra +
                                         w_cov * w_cov;

                    const StateReal dV = m.sqrtG * (R2 * dxi * dxi) * dz;

                    mass    += (long double)rho * (long double)dV;
                    theta_int += (long double)rho_theta * (long double)dV;
                    ke      += (long double)(0.5 * rho * v2) * (long double)dV;
                    ie      += (long double)(rho * gas_constants::c_v * T) *
                               (long double)dV;
                    pe      += (long double)(rho * Phi) * (long double)dV;

                    // 角动量：ρ (u r cosφ + Ω r² cos²φ)
                    // 用经度 ξ 方向速度 u（反变经向）与纬度 cosφ
                    const StateReal coslat = std::cos(m.lat);
                    const StateReal r_local = radius + z;
                    const StateReal u_zon = u_contra;  // ξ 方向近似纬向
                    angmom += (long double)(rho * (u_zon * r_local * coslat +
                        omega * r_local * r_local * coslat * coslat)) *
                        (long double)dV;
                }
            }
        }
    }

    snap.total_mass = static_cast<StateReal>(mass);
    snap.total_theta = static_cast<StateReal>(theta_int);
    snap.kinetic_energy = static_cast<StateReal>(ke);
    snap.internal_energy = static_cast<StateReal>(ie);
    snap.potential_energy = static_cast<StateReal>(pe);
    snap.total_energy = snap.kinetic_energy + snap.internal_energy +
                        snap.potential_energy;
    snap.total_angular_momentum = static_cast<StateReal>(angmom);
    return snap;
}

std::vector<StateReal> ConservationDiagnostics::relative_drift(
    const ConservationSnapshot& initial,
    const ConservationSnapshot& current) {
    std::vector<StateReal> drift(4, 0.0);
    auto rel = [](StateReal c, StateReal i) {
        if (std::abs(i) < 1e-30) return 0.0;
        return std::abs(c - i) / std::abs(i);
    };
    drift[0] = rel(current.total_mass, initial.total_mass);
    drift[1] = rel(current.total_theta, initial.total_theta);
    drift[2] = rel(current.total_energy, initial.total_energy);
    drift[3] = rel(current.total_angular_momentum, initial.total_angular_momentum);
    return drift;
}

}  // namespace cubed_sph::diagnostics
