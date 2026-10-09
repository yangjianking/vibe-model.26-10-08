// =============================================================================
//  状态实现
// =============================================================================

#include "cubed_sph/dynamics/state.hpp"

namespace cubed_sph::dynamics {

using grid::kNumPanels;

State::State(const grid::CubedSphereGrid& grid, int nlev)
    : grid_(grid), nlev_(nlev) {
    ncell2d_ = static_cast<IIndex>(kNumPanels) * grid.panel_n() * grid.panel_n();
    const IIndex total = ncell2d_ * nlev_;
    rho.assign(total, 0.0);
    rho_u.assign(total, 0.0);
    rho_v.assign(total, 0.0);
    rho_w.assign(total, 0.0);
    rho_theta.assign(total, 0.0);
    exner.assign(total, 0.0);
    theta.assign(total, 0.0);
    geopotential.assign(total, 0.0);
}

State& State::operator=(const State& other) {
    // 引用成员 grid_ 保持不变（两侧应绑定同一网格）。
    // nlev_ 与 ncell2d_ 亦不覆盖，仅拷贝字段数据。
    rho          = other.rho;
    rho_u        = other.rho_u;
    rho_v        = other.rho_v;
    rho_w        = other.rho_w;
    rho_theta    = other.rho_theta;
    exner        = other.exner;
    theta        = other.theta;
    geopotential = other.geopotential;
    return *this;
}

IIndex State::index(int panel, int i, int j, int k) const {
    const int n = grid_.panel_n();
    return static_cast<IIndex>(panel) * n * n * nlev_ +
           static_cast<IIndex>(k) * n * n +
           static_cast<IIndex>(j) * n + i;
}

void State::diagnose_velocities() {
    // 由动量密度还原速度场（供诊断与 I/O）。完整实现见 diagnostics.cpp。
    // u = ρu / ρ（协变分量）
}

}  // namespace cubed_sph::dynamics
