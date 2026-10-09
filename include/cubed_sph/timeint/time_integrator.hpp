// =============================================================================
//  时间积分（Time Integration）—— 半隐式时间推进
// =============================================================================
//  采用时间分裂（time splitting，提示词 P2 第 3 条）：
//    - 慢过程（平流）：守恒型有限体积平流，大时间步 Δt_adv 仅受精度约束；
//    - 快波（重力波、声波）：半隐式（SI）——对产生快波的线性项用
//      Crank–Nicolson 隐式平均（off-centering 参数 α ∈ [0.5, 0.6]），
//      其余项显式。
//
//  半隐式推导（从线性化方程到 Helmholtz 方程）：
//    线性化快波项（气压梯度力与散度耦合、浮力与垂直速度耦合）用
//    Crank–Nicolson 时间平均后，可消元得到关于 Exner 压力扰动 δπ' 的
//    三维椭圆方程：
//      (I - β² Δt² c_s² ∇²) δπ' = RHS
//    其中 β = α Δt 为隐式权重相关的系数，c_s 为参考态声速。
//
//  完整推导（连续→离散、截断误差、von Neumann 稳定性）见
//    docs/numerics/semi_implicit.md 与 docs/theory/helmholtz.md
//
//  参考文献：
//    Robert (1981), Atmos.-Ocean 19, 35–52（半隐式方案）
//    Tanguay et al. (1990), Mon. Wea. Rev. 118, 1970–1980（可压缩半隐式）
//    Staniforth & Côté (1991), Mon. Wea. Rev. 119, 2206–2223（综述）
//    Wood et al. (2014), QJRMS 140, 1505–1528（ENDGame 半隐式）
// =============================================================================

#pragma once

#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/dynamics/equations.hpp"
#include "cubed_sph/grid/vertical_coordinate.hpp"
#include "cubed_sph/common/config.hpp"

namespace cubed_sph::timeint {

// ---------------------------------------------------------------------------
// 时间积分参数（由配置注入）
// ---------------------------------------------------------------------------
struct TimeStepperParams {
    StateReal dt = 300.0;         // 大时间步 [s]
    StateReal off_centering = 0.5; // α ∈ [0.5, 0.6]
    int n_split = 1;              // 声波子步数（显式子步，可选）
    bool use_semi_implicit = true;
};

// ---------------------------------------------------------------------------
// 时间积分器抽象接口（可扩展：半隐式 Euler、IMEX-RK/HEVI 等）
// ---------------------------------------------------------------------------
class TimeIntegrator {
public:
    virtual ~TimeIntegrator() = default;

    // 推进一个时间步：state 原地更新
    virtual void step(dynamics::State& state) = 0;

    // 返回方案名与当前步数（用于日志/provenance）
    virtual const char* name() const = 0;
    virtual IIndex step_count() const = 0;
};

// ---------------------------------------------------------------------------
// 半隐式时间积分器
// ---------------------------------------------------------------------------
// 完整时间分裂推进（提示词 P2 第 3 条）：
//   1. 慢过程（平流）显式：大时间步 Δt_adv，仅受精度约束；
//   2. 快波（重力波/声波）半隐式：Crank–Nicolson 隐式平均，求解 Helmholtz
//      方程 (I - β²Δt²c_s²∇²) δπ' = RHS；
//   3. 物理倾向 operator splitting 叠加。
class SemiImplicitIntegrator : public TimeIntegrator {
public:
    SemiImplicitIntegrator(dynamics::DynamicsCore& core,
                           const grid::VerticalCoordinate& vert,
                           const TimeStepperParams& params);

    void step(dynamics::State& state) override;
    const char* name() const override { return "SemiImplicit"; }
    IIndex step_count() const override { return nsteps_; }

private:
    dynamics::DynamicsCore& core_;
    const grid::VerticalCoordinate& vert_;
    TimeStepperParams params_;
    IIndex nsteps_ = 0;

    // 求解半隐式 Helmholtz 方程（隐式快波修正），返回 δπ' 迭代情况
    int solve_helmholtz(dynamics::State& state);
};

// ---------------------------------------------------------------------------
// 显式时间积分器（对照/调试用，SSP-RK3 三阶段）
// ---------------------------------------------------------------------------
class ExplicitRK3Integrator : public TimeIntegrator {
public:
    ExplicitRK3Integrator(dynamics::DynamicsCore& core,
                          const TimeStepperParams& params);

    void step(dynamics::State& state) override;
    const char* name() const override { return "ExplicitRK3"; }
    IIndex step_count() const override { return nsteps_; }

private:
    dynamics::DynamicsCore& core_;
    TimeStepperParams params_;
    IIndex nsteps_ = 0;

    // SSP-RK3 单个阶段：q ← q + c * Δt * F(q)（低存储形式）
    void rk3_stage(dynamics::State& state, StateReal c, dynamics::State& work);
};

}  // namespace cubed_sph::timeint
