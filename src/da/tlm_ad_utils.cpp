// =============================================================================
//  状态空间线性代数工具（内积 / axpy / 随机扰动）
// =============================================================================
//  从 tlm_ad.cpp 中抽离的纯工具函数，不依赖 DynamicsCore/Config，
//  供 da 模块内部（TL/AD、4D-Var）与单元测试独立链接复用。
// =============================================================================

#include "cubed_sph/da/tlm_ad.hpp"

#include <algorithm>

namespace cubed_sph::da {

StateReal state_inner_product(const dynamics::State& a,
                              const dynamics::State& b) {
    // 长双精度补偿求和，避免大状态空间下浮点求和误差影响点积检验精度。
    const IIndex n = a.ncell();
    long double acc = 0.0L;
    for (IIndex i = 0; i < n; ++i) {
        acc += (long double)a.rho[i]        * (long double)b.rho[i];
        acc += (long double)a.rho_u[i]      * (long double)b.rho_u[i];
        acc += (long double)a.rho_v[i]      * (long double)b.rho_v[i];
        acc += (long double)a.rho_w[i]      * (long double)b.rho_w[i];
        acc += (long double)a.rho_theta[i]  * (long double)b.rho_theta[i];
        acc += (long double)a.exner[i]      * (long double)b.exner[i];
        acc += (long double)a.theta[i]      * (long double)b.theta[i];
        acc += (long double)a.geopotential[i] * (long double)b.geopotential[i];
    }
    return static_cast<StateReal>(acc);
}

void state_axpy(StateReal alpha, const dynamics::State& x, dynamics::State& y) {
    const IIndex n = x.ncell();
    for (IIndex i = 0; i < n; ++i) {
        y.rho[i]        += alpha * x.rho[i];
        y.rho_u[i]      += alpha * x.rho_u[i];
        y.rho_v[i]      += alpha * x.rho_v[i];
        y.rho_w[i]      += alpha * x.rho_w[i];
        y.rho_theta[i]  += alpha * x.rho_theta[i];
        y.exner[i]      += alpha * x.exner[i];
        y.theta[i]      += alpha * x.theta[i];
        y.geopotential[i] += alpha * x.geopotential[i];
    }
}

void state_random(dynamics::State& s, std::uint64_t seed, StateReal scale) {
    // xorshift64 可复现随机数
    auto next = [&seed]() -> StateReal {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        const std::uint64_t v = seed;
        return (static_cast<StateReal>(v >> 11) / 9007199254740992.0) - 1.0;
    };
    const IIndex n = s.ncell();
    for (IIndex i = 0; i < n; ++i) {
        s.rho[i]        = scale * next();
        s.rho_u[i]      = scale * next();
        s.rho_v[i]      = scale * next();
        s.rho_w[i]      = scale * next();
        s.rho_theta[i]  = scale * next();
        s.exner[i]      = scale * next();
        s.theta[i]      = scale * next();
        s.geopotential[i] = scale * next();
    }
}

}  // namespace cubed_sph::da
