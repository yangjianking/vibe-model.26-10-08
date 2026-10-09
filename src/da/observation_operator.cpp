// =============================================================================
//  观测算子与误差协方差实现
// =============================================================================

#include "cubed_sph/da/observation_operator.hpp"

#include <cmath>
#include <algorithm>

namespace cubed_sph::da {

using grid::kNumPanels;

// =============================================================================
// DiagonalBackgroundCovariance
// =============================================================================
StateReal DiagonalBackgroundCovariance::var_for_field(int field_id) const {
    // field_id: 0=rho, 1..3=动量(rho_u/rho_v/rho_w), 4=rho_theta
    //           5..7=诊断(exner/theta/geopotential)
    switch (field_id) {
        case 0:  return sigma_rho * sigma_rho;
        case 1: case 2: case 3: return sigma_mom * sigma_mom;
        case 4:  return sigma_theta * sigma_theta;
        default: return 0.0;  // 诊断量不作背景约束
    }
}

void DiagonalBackgroundCovariance::apply_B(const dynamics::State& v,
                                           dynamics::State& Bv) const {
    const IIndex n = v.ncell();
    const StateReal vr = var_for_field(0);
    const StateReal vm = var_for_field(1);
    const StateReal vt = var_for_field(4);
    for (IIndex i = 0; i < n; ++i) {
        Bv.rho[i]        = vr * v.rho[i];
        Bv.rho_u[i]      = vm * v.rho_u[i];
        Bv.rho_v[i]      = vm * v.rho_v[i];
        Bv.rho_w[i]      = vm * v.rho_w[i];
        Bv.rho_theta[i]  = vt * v.rho_theta[i];
        Bv.exner[i]      = 0.0;
        Bv.theta[i]      = 0.0;
        Bv.geopotential[i] = 0.0;
    }
}

void DiagonalBackgroundCovariance::apply_Binv(const dynamics::State& v,
                                              dynamics::State& Binvv) const {
    const IIndex n = v.ncell();
    const StateReal ir = 1.0 / var_for_field(0);
    const StateReal im = 1.0 / var_for_field(1);
    const StateReal it = 1.0 / var_for_field(4);
    for (IIndex i = 0; i < n; ++i) {
        Binvv.rho[i]        = ir * v.rho[i];
        Binvv.rho_u[i]      = im * v.rho_u[i];
        Binvv.rho_v[i]      = im * v.rho_v[i];
        Binvv.rho_w[i]      = im * v.rho_w[i];
        Binvv.rho_theta[i]  = it * v.rho_theta[i];
        Binvv.exner[i]      = 0.0;
        Binvv.theta[i]      = 0.0;
        Binvv.geopotential[i] = 0.0;
    }
}

void DiagonalBackgroundCovariance::sqrt_B(const dynamics::State& w,
                                          dynamics::State& B12w) const {
    const IIndex n = w.ncell();
    const StateReal sr = std::sqrt(var_for_field(0));
    const StateReal sm = std::sqrt(var_for_field(1));
    const StateReal st = std::sqrt(var_for_field(4));
    for (IIndex i = 0; i < n; ++i) {
        B12w.rho[i]        = sr * w.rho[i];
        B12w.rho_u[i]      = sm * w.rho_u[i];
        B12w.rho_v[i]      = sm * w.rho_v[i];
        B12w.rho_w[i]      = sm * w.rho_w[i];
        B12w.rho_theta[i]  = st * w.rho_theta[i];
        B12w.exner[i]      = 0.0;
        B12w.theta[i]      = 0.0;
        B12w.geopotential[i] = 0.0;
    }
}

void DiagonalBackgroundCovariance::sqrt_Binv(const dynamics::State& v,
                                             dynamics::State& Bmin12v) const {
    const IIndex n = v.ncell();
    const StateReal ir = 1.0 / std::sqrt(var_for_field(0));
    const StateReal im = 1.0 / std::sqrt(var_for_field(1));
    const StateReal it = 1.0 / std::sqrt(var_for_field(4));
    for (IIndex i = 0; i < n; ++i) {
        Bmin12v.rho[i]        = ir * v.rho[i];
        Bmin12v.rho_u[i]      = im * v.rho_u[i];
        Bmin12v.rho_v[i]      = im * v.rho_v[i];
        Bmin12v.rho_w[i]      = im * v.rho_w[i];
        Bmin12v.rho_theta[i]  = it * v.rho_theta[i];
        Bmin12v.exner[i]      = 0.0;
        Bmin12v.theta[i]      = 0.0;
        Bmin12v.geopotential[i] = 0.0;
    }
}

// =============================================================================
// DiagonalObservationCovariance
// =============================================================================
DiagonalObservationCovariance::DiagonalObservationCovariance(
    const std::vector<io::Observation>& obs) {
    var_.reserve(obs.size());
    for (const auto& o : obs) {
        var_.push_back(o.error * o.error);
    }
}

void DiagonalObservationCovariance::apply_R(const std::vector<StateReal>& v,
                                            std::vector<StateReal>& Rv) const {
    Rv.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        Rv[i] = var_[i] * v[i];
    }
}

void DiagonalObservationCovariance::apply_Rinv(const std::vector<StateReal>& v,
                                               std::vector<StateReal>& Rinvv) const {
    Rinvv.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        Rinvv[i] = v[i] / var_[i];
    }
}

// =============================================================================
// SimpleThetaObservationOperator
// =============================================================================
SimpleThetaObservationOperator::SimpleThetaObservationOperator(
    const std::vector<io::Observation>& obs,
    const grid::CubedSphereGrid& grid,
    int nlev)
    : obs_(obs), grid_(grid), nlev_(nlev), nobs_(obs.size()) {}

void SimpleThetaObservationOperator::forward(const dynamics::State& state,
                                             std::vector<StateReal>& hx) const {
    // 占位：假设观测位置最近邻映射到某网格点，取该点 θ。
    // 完整实现需水平/垂直插值。此处给出可编译、满足线性/伴随一致性的骨架：
    // 观测 i 映射到固定网格单元 index = (i mod ncell2d)，取 state.theta[index]。
    hx.assign(nobs_, 0.0);
    const int n = grid_.panel_n();
    const IIndex ncell2d = static_cast<IIndex>(kNumPanels) * n * n;
    for (size_t i = 0; i < nobs_; ++i) {
        const IIndex flat = static_cast<IIndex>(i) % ncell2d;
        hx[i] = state.theta[flat];  // 简化：取表层（k=0）theta
    }
}

void SimpleThetaObservationOperator::tangent_linear(
    const dynamics::State& dx, std::vector<StateReal>& hdx) const {
    hdx.assign(nobs_, 0.0);
    const int n = grid_.panel_n();
    const IIndex ncell2d = static_cast<IIndex>(kNumPanels) * n * n;
    for (size_t i = 0; i < nobs_; ++i) {
        const IIndex flat = static_cast<IIndex>(i) % ncell2d;
        hdx[i] = dx.theta[flat];
    }
}

void SimpleThetaObservationOperator::adjoint(const std::vector<StateReal>& dy,
                                             dynamics::State& atdx) const {
    const int n = grid_.panel_n();
    const IIndex ncell2d = static_cast<IIndex>(kNumPanels) * n * n;
    for (size_t i = 0; i < dy.size() && i < nobs_; ++i) {
        const IIndex flat = static_cast<IIndex>(i) % ncell2d;
        atdx.theta[flat] += dy[i];  // 累加（伴随为转置 + 散射）
    }
}

}  // namespace cubed_sph::da
