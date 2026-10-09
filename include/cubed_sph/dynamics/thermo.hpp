// =============================================================================
//  热力学状态（Thermodynamics）—— 理想气体关系与 Exner 压力
// =============================================================================
//  干空气理想气体状态方程与热力学关系：
//    p = ρ R_d T                                  （理想气体）
//    θ = T (p0/p)^κ                               （位温定义）
//    π = (p/p0)^κ = R_d θ / (R_d θ) ...          （Exner 压力）
//    Φ = g z                                      （重力位）
//
//  其中 R_d = 287.05 J/(kg·K)（干空气气体常数），
//        c_p = 1004.5 J/(kg·K)（定压比热），
//        κ  = R_d / c_p ≈ 0.2857，
//        p0 = 100000 Pa（参考气压），
//        g  = 9.80665 m/s²（重力加速度）。
//
//  预后变量若采用 ρθ（位温密度），则温度、气压、Exner 压力均由诊断关系
//  恢复（Harris et al. 2021 §4）。
//
//  参考文献：
//    Harris et al. (2021), GFDL TM GFDL2021001 §4（热力学关系）
//    Dutton, J. A. (1986), The Ceaseless Wind（大气热力学基础）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"

namespace cubed_sph::dynamics {

// 干空气物理常数（国际标准，单位 SI）
namespace gas_constants {
    inline constexpr StateReal R_d  = 287.05;        // 干空气气体常数 [J/(kg·K)]
    inline constexpr StateReal c_p  = 1004.5;        // 定压比热 [J/(kg·K)]
    inline constexpr StateReal c_v  = 717.5;         // 定容比热 [J/(kg·K)]
    inline constexpr StateReal kappa = R_d / c_p;    // ≈ 0.2857
    inline constexpr StateReal p0   = 100000.0;      // 参考气压 [Pa]
    inline constexpr StateReal g    = 9.80665;       // 重力加速度 [m/s²]
}

// ---------------------------------------------------------------------------
// 位温 → 温度（由 Exner 压力）： T = θ π
// ---------------------------------------------------------------------------
inline StateReal theta_to_temperature(StateReal theta, StateReal exner) {
    return theta * exner;
}

// ---------------------------------------------------------------------------
// 密度与位温密度 → Exner 压力
//    π = (ρ θ R_d / p0)^{κ/(1-κ)}   （由 ρθ 与 ρ 恢复 π）
// 推导见 docs/theory/thermodynamics.md；此处用标准闭合关系。
// ---------------------------------------------------------------------------
inline StateReal exner_from_rho_theta(StateReal rho, StateReal rho_theta) {
    const StateReal theta = rho_theta / rho;
    const StateReal base = theta * gas_constants::R_d / gas_constants::p0;
    return std::pow(base, gas_constants::kappa / (1.0 - gas_constants::kappa));
}

// ---------------------------------------------------------------------------
// Exner 压力 → 气压： p = p0 π^{1/κ}
// ---------------------------------------------------------------------------
inline StateReal exner_to_pressure(StateReal exner) {
    return gas_constants::p0 * std::pow(exner, 1.0 / gas_constants::kappa);
}

// ---------------------------------------------------------------------------
// 声速（用于 Helmholtz 算子与 CFL 约束）
//    c_s² = γ R_d T = c_p/c_v · R_d T
// ---------------------------------------------------------------------------
inline StateReal sound_speed_sq(StateReal temperature) {
    const StateReal gamma = gas_constants::c_p / gas_constants::c_v;
    return gamma * gas_constants::R_d * temperature;
}

}  // namespace cubed_sph::dynamics
