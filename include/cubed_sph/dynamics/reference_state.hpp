// =============================================================================
//  参考态（Reference State）—— 静力平衡剖面
// =============================================================================
//  半隐式时间推进需要将动力方程在静力平衡参考态附近线性化。参考态满足
//  静力平衡：
//
//    ∂p̄/∂z = -ρ̄ g          （静力平衡）
//    p̄ = ρ̄ R_d T̄           （理想气体）
//
//  采用等温参考态（isothermal reference atmosphere，T̄ = const），解析解为：
//
//    p̄(z) = p_s exp(-z / H),   ρ̄(z) = ρ_s exp(-z / H),   H = R_d T̄ / g
//
//  其中 H 为标高（scale height），p_s、ρ_s 为地面参考值。等温参考态下
//  Exner 压力、声速均有解析表达，极大简化半隐式 Helmholtz 算子的构造。
//
//  预后变量的扰动分解（Harris et al. 2021 §4）：
//    ρ  = ρ̄ + ρ',   ρθ = ρ̄θ̄ + (ρθ)',   π = π̄ + π'
//  快波线性项（气压梯度 ↔ 散度、浮力 ↔ 垂直速度）用参考态系数展开，
//  半隐式隐式部分消元后得到关于 π' 的 Helmholtz 方程。
//
//  参考文献：
//    Harris et al. (2021), GFDL TM GFDL2021001, §4–5
//    Wood et al. (2014), QJRMS 140, 1505–1528（ENDGame 参考态）
//    Staniforth & Côté (1991), Mon. Wea. Rev. 119, 2206–2223（半隐式线性化）
// =============================================================================

#pragma once

#include "cubed_sph/grid/vertical_coordinate.hpp"
#include "cubed_sph/dynamics/thermo.hpp"
#include "cubed_sph/common/types.hpp"

#include <vector>

namespace cubed_sph::dynamics {

// ---------------------------------------------------------------------------
// 参考态（等温静力平衡）
// ---------------------------------------------------------------------------
class ReferenceState {
public:
    // 用等温参考温度 T_ref 与地面参考气压 p_s 构造参考态
    ReferenceState(const grid::VerticalCoordinate& vert,
                   StateReal T_ref = 250.0,
                   StateReal p_surf = 100000.0);

    // ---- 参考态剖面（按垂直主层） ----
    StateReal pressure(int k) const { return p_[k]; }       // p̄ [Pa]
    StateReal density(int k) const { return rho_[k]; }      // ρ̄ [kg/m³]
    StateReal temperature() const { return T_ref_; }        // T̄ [K]（等温）
    StateReal exner(int k) const { return exner_[k]; }      // π̄
    StateReal theta(int k) const { return theta_[k]; }      // θ̄ [K]
    StateReal sound_speed_sq(int k) const { return cs2_[k]; } // c_s² [m²/s²]
    StateReal height(int k) const { return z_[k]; }          // 主层高度 [m]

    // ---- 派生量 ----
    StateReal scale_height() const { return H_; }  // 标高 [m]

private:
    const grid::VerticalCoordinate& vert_;
    StateReal T_ref_;
    StateReal p_surf_;
    StateReal H_;  // 标高 = R_d T̄ / g

    std::vector<StateReal> z_;      // 主层高度
    std::vector<StateReal> p_;      // 参考气压
    std::vector<StateReal> rho_;    // 参考密度
    std::vector<StateReal> exner_;  // 参考 Exner 压力
    std::vector<StateReal> theta_;  // 参考位温
    std::vector<StateReal> cs2_;    // 参考声速平方
};

}  // namespace cubed_sph::dynamics
