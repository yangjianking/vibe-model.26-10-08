// =============================================================================
//  参考态实现（等温静力平衡剖面）
// =============================================================================

#include "cubed_sph/dynamics/reference_state.hpp"

#include <cmath>

namespace cubed_sph::dynamics {

using gas_constants::R_d;
using gas_constants::c_p;
using gas_constants::c_v;
using gas_constants::g;
using gas_constants::p0;
using gas_constants::kappa;

ReferenceState::ReferenceState(const grid::VerticalCoordinate& vert,
                               StateReal T_ref, StateReal p_surf)
    : vert_(vert), T_ref_(T_ref), p_surf_(p_surf) {
    // 标高 H = R_d T̄ / g
    H_ = R_d * T_ref_ / g;

    const int nlev = vert.nlev();
    z_.resize(nlev);
    p_.resize(nlev);
    rho_.resize(nlev);
    exner_.resize(nlev);
    theta_.resize(nlev);
    cs2_.resize(nlev);

    // 地面参考密度（由状态方程）
    const StateReal rho_surf = p_surf_ / (R_d * T_ref_);
    // 声速平方（等温）：c_s² = γ R_d T̄
    const StateReal gamma = c_p / c_v;
    const StateReal cs2 = gamma * R_d * T_ref_;

    for (int k = 0; k < nlev; ++k) {
        z_[k] = vert.z_full(k);
        // 等温静力平衡：p̄(z) = p_s exp(-z/H)
        p_[k]   = p_surf_ * std::exp(-z_[k] / H_);
        rho_[k] = rho_surf * std::exp(-z_[k] / H_);
        // Exner 压力 π̄ = (p̄/p0)^κ
        exner_[k] = std::pow(p_[k] / p0, kappa);
        // 参考位温 θ̄ = T̄ / π̄
        theta_[k] = T_ref_ / exner_[k];
        cs2_[k] = cs2;
    }
}

}  // namespace cubed_sph::dynamics
