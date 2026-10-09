// =============================================================================
//  平流模块（Advection）—— 通量形式平流
// =============================================================================
//  水平平流采用守恒型有限体积 + 高阶重构（PPM / 三阶上游），保证正性
//  （qc、qv ≥ 0）与守恒。本阶段提供平流算子的接口与守恒框架，具体重构
//  方案（Lin–Rood 1996 通量形式半拉格朗日 vs Wicker–Skamarock 2002 三阶）
//  作为可扩展策略注入。
//
//  通量形式： ∂(ρφ)/∂t = -∇·(ρφ v)，其中 φ 为示踪量或单位质量量。
//
//  参考文献：
//    Lin & Rood (1996), Mon. Wea. Rev. 124, 2046–2070（通量形式半拉格朗日）
//    Wicker & Skamarock (2002), Mon. Wea. Rev. 130, 2088–2097（三阶上游）
//    Colella & Woodward (1984), J. Comput. Phys. 54, 174–201（PPM）
// =============================================================================

#pragma once

#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/grid/differential_operators.hpp"

#include <vector>

namespace cubed_sph::dynamics {

// ---------------------------------------------------------------------------
// 平流方案选择（编译期/配置注入的策略）
// ---------------------------------------------------------------------------
enum class AdvectionScheme {
    Center2,     // 二阶中心差分（骨架/调试）
    Upwind3,     // 三阶上游（Wicker–Skamarock 2002）
    PPM,         // 分段抛物（Colella & Woodward 1984）
};

// ---------------------------------------------------------------------------
// 平流算子接口（可扩展：新方案实现 compute_flux 即可）
// ---------------------------------------------------------------------------
class Advection {
public:
    Advection(AdvectionScheme scheme, const grid::CubedSphereGrid& grid);

    // 计算质量/示踪量平流倾向： -∇·(ρ φ v)
    // 输入：密度场 rho、被平流量 phi、协变速度 u_cov/v_cov
    // 输出：倾向 d(ρφ)/dt 累加到 tendency
    void advect(const State& state,
                const std::vector<StateReal>& phi,   // 被平流标量
                const std::vector<StateReal>& u_cov,
                const std::vector<StateReal>& v_cov,
                std::vector<StateReal>& tendency);

    // 返回当前方案名（用于日志与 provenance）
    const char* scheme_name() const;

private:
    AdvectionScheme scheme_;
    const grid::CubedSphereGrid& grid_;
};

}  // namespace cubed_sph::dynamics
