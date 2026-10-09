// =============================================================================
//  垂直坐标实现
// =============================================================================
//  构建垂直层级与度量项。地形跟随坐标：
//    z(η) = z_top * η + h_surf * b(η)
//  其中 b(η) 为 Schär 2002 平滑地形衰减函数。
// =============================================================================

#include "cubed_sph/grid/vertical_coordinate.hpp"

#include <cmath>

namespace cubed_sph::grid {

// ---------------------------------------------------------------------------
// Schär 2002 地形衰减函数：
//   b(η) = sinh{(η_top - η) / s} / sinh{η_top / s}
// 满足 b(0)=1（地面），b(η_top)=0（模式顶）。s 为衰减尺度（越小地形衰减越快，
// 地形影响越局限于近地面）。
// ---------------------------------------------------------------------------
StateReal VerticalCoordinate::decay(StateReal eta, StateReal s) {
    // η_top = 1（归一化坐标的模式顶）
    const StateReal eta_top = 1.0;
    return std::sinh((eta_top - eta) / s) / std::sinh(eta_top / s);
}

StateReal VerticalCoordinate::decay_deriv(StateReal eta, StateReal s) {
    const StateReal eta_top = 1.0;
    // d b/dη = -(1/s) cosh{(η_top-η)/s} / sinh{η_top/s}
    return -std::cosh((eta_top - eta) / s) / (s * std::sinh(eta_top / s));
}

VerticalCoordinate::VerticalCoordinate(const VerticalParams& params)
    : params_(params),
      nlev_(params.nlev),
      z_top_(params.z_top),
      type_(params.type),
      stagger_(params.stagger) {
    z_full_.assign(nlev_, 0.0);
    z_half_.assign(nlev_ + 1, 0.0);
    jacobian_.assign(nlev_, 0.0);
    deta_dz_.assign(nlev_ + 1, 0.0);
    build_levels();
}

void VerticalCoordinate::rebuild() {
    build_levels();
}

void VerticalCoordinate::build_levels() {
    // 归一化 η 坐标：主层中心 η_k，半层界面 η_{k+1/2}
    // 半层等距：η_{k+1/2} = k / nlev,  k = 0..nlev
    // 主层中心：η_k = (k + 0.5) / nlev,  k = 0..nlev-1
    const StateReal d_eta = 1.0 / static_cast<StateReal>(nlev_);

    for (int k = 0; k <= nlev_; ++k) {
        const StateReal eta_half = static_cast<StateReal>(k) * d_eta;

        if (type_ == VerticalCoordType::TerrainFollowing) {
            // 地形跟随：z = z_top * η + h_surf * b(η)
            z_half_[k] = z_top_ * eta_half + h_surf_ * decay(eta_half, params_.s);
        } else {
            // 纯高度：z = z_top * η
            z_half_[k] = z_top_ * eta_half;
        }
    }

    for (int k = 0; k < nlev_; ++k) {
        const StateReal eta_full = (static_cast<StateReal>(k) + 0.5) * d_eta;

        if (type_ == VerticalCoordType::TerrainFollowing) {
            z_full_[k] = z_top_ * eta_full + h_surf_ * decay(eta_full, params_.s);
            // J = ∂z/∂η = z_top + h_surf * db/dη
            jacobian_[k] = z_top_ + h_surf_ * decay_deriv(eta_full, params_.s);
        } else {
            z_full_[k] = z_top_ * eta_full;
            jacobian_[k] = z_top_;  // 纯高度下 ∂z/∂η = z_top
        }
    }

    // 半层 ∂η/∂z = 1/J（用相邻主层 Jacobian 的调和/算术平均，此处取算术平均）
    for (int k = 1; k < nlev_; ++k) {
        deta_dz_[k] = 1.0 / (0.5 * (jacobian_[k - 1] + jacobian_[k]));
    }
    deta_dz_[0] = 1.0 / jacobian_[0];
    deta_dz_[nlev_] = 1.0 / jacobian_[nlev_ - 1];
}

}  // namespace cubed_sph::grid
