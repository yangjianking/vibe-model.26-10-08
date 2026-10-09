// =============================================================================
//  切线性模型（TL）与伴随模型（AD）（P6 手写 TL/AD）
// =============================================================================
//  4D-Var 需要预报模型 M 的切线性算子 M' 与伴随算子 Mᵀ：
//    - 切线性（Tangent Linear, TL）：δx(t₀) → δx(t_k)，即 M_k δx
//    - 伴随（Adjoint, AD）：δx̂(t_k) → δx̂(t₀)，即 M_kᵀ δx̂
//
//  二者满足点积检验（dot-product test，用于验证互为转置）：
//      ⟨M'u, v⟩ = ⟨u, Mᵀv⟩  ∀ u, v
//  数值上要求相对误差 ≲ 1e-12（fp64）。
//
//  切线性模型是预报模型 M 对状态的 Gateaux 导数（线性化）：
//      M'(x) = lim_{ε→0} [M(x + ε δx) - M(x)] / ε
//  伴随模型是 M' 的转置（对欧氏内积），用于反向传播梯度。
//
//  本模块定义 TLM 与 ADM 的抽象接口，并提供一个基于有限差分的
//  "动态 TLM"（用于点积检验的参照，非生产路径）与一个手写线性化骨架。
//
//  参考文献：
//    Errico (1997), Bull. Amer. Meteor. Soc. 78, 2577–2591（TL/AD 综述）
//    Giering & Kaminski (1998), ACM TOMS 24, 437–474（自动微分）
//    Lorenc (2003), QJRMS 129, 3183–3203（增量 4D-Var 的 TL/AD）
// =============================================================================

#pragma once

#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/dynamics/equations.hpp"
#include "cubed_sph/grid/vertical_coordinate.hpp"
#include "cubed_sph/timeint/time_integrator.hpp"
#include "cubed_sph/common/config.hpp"

namespace cubed_sph::da {

// ---------------------------------------------------------------------------
// 切线性模型接口：给定轨迹 x_ref，计算 M'(x_ref) δx
// ---------------------------------------------------------------------------
class TangentLinearModel {
public:
    virtual ~TangentLinearModel() = default;

    // 沿参考轨迹 x_ref，将扰动 dx 从 t₀ 传播到 t₁（dx 原地更新）
    virtual void propagate(const dynamics::State& x_ref,
                           dynamics::State& dx,
                           StateReal t0, StateReal t1) = 0;

    virtual const char* name() const = 0;
};

// ---------------------------------------------------------------------------
// 伴随模型接口：给定轨迹 x_ref，反向传播 dx̂（从 t₁ 到 t₀）
// ---------------------------------------------------------------------------
class AdjointModel {
public:
    virtual ~AdjointModel() = default;

    // 沿参考轨迹 x_ref，将伴随量 dx̂ 从 t₁ 反传到 t₀（dx̂ 原地更新）
    virtual void propagate(const dynamics::State& x_ref,
                           dynamics::State& dx_hat,
                           StateReal t0, StateReal t1) = 0;

    virtual const char* name() const = 0;
};

// ---------------------------------------------------------------------------
// 点积检验（dot-product test）
// ---------------------------------------------------------------------------
// 验证 TL 与 AD 互为转置。给定随机扰动 u（状态空间）、v（状态空间），检查
//      ⟨M'u, v⟩ ≈ ⟨u, Mᵀv⟩
// 返回相对误差。要求 ≲ 1e-12（fp64）。
// ---------------------------------------------------------------------------
class DotProductTest {
public:
    DotProductTest(TangentLinearModel& tlm, AdjointModel& adm);

    // 执行一次点积检验，返回相对误差与左右两侧点积值（供诊断）
    StateReal run(const dynamics::State& x_ref,
                  StateReal t0, StateReal t1,
                  StateReal& lhs, StateReal& rhs);

    // 是否通过（默认容差 1e-10）
    static bool passed(StateReal rel_err, StateReal tol = 1e-10);

private:
    TangentLinearModel& tlm_;
    AdjointModel& adm_;
};

// ---------------------------------------------------------------------------
// 有限差分切线性模型（点积检验参照，非生产路径）
// ---------------------------------------------------------------------------
// 用中心差分近似 M'：
//      M'(x) δx ≈ [M(x + ε δx) - M(x - ε δx)] / (2ε)
// 仅用于验证手写 TL 的正确性，不可用于生产（昂贵且受舍入误差影响）。
// ---------------------------------------------------------------------------
class FiniteDifferenceTLM : public TangentLinearModel {
public:
    FiniteDifferenceTLM(dynamics::DynamicsCore& core,
                        const grid::VerticalCoordinate& vert,
                        StateReal dt, StateReal eps);

    void propagate(const dynamics::State& x_ref, dynamics::State& dx,
                   StateReal t0, StateReal t1) override;
    const char* name() const override { return "FiniteDifferenceTLM"; }

private:
    dynamics::DynamicsCore& core_;
    const grid::VerticalCoordinate& vert_;
    StateReal dt_;
    StateReal eps_;
};

// ---------------------------------------------------------------------------
// 手写切线性模型（单时间步、线性化动力倾向）
// ---------------------------------------------------------------------------
// 对显式动力倾向 F(x) 线性化：δx ← δx + Δt · F'(x) δx
// F' 为 F 的 Jacobian，手写实现（此处给出占位：线性化平流 + 线性化气压梯度
// 的骨架，完整推导见 docs/theory/tlm_ad.md）。
// ---------------------------------------------------------------------------
class HandwrittenTLM : public TangentLinearModel {
public:
    HandwrittenTLM(dynamics::DynamicsCore& core,
                   const grid::VerticalCoordinate& vert,
                   StateReal dt);

    void propagate(const dynamics::State& x_ref, dynamics::State& dx,
                   StateReal t0, StateReal t1) override;
    const char* name() const override { return "HandwrittenTLM"; }

private:
    dynamics::DynamicsCore& core_;
    const grid::VerticalCoordinate& vert_;
    StateReal dt_;
};

// ---------------------------------------------------------------------------
// 手写伴随模型（单时间步、转置线性化动力倾向）
// ---------------------------------------------------------------------------
// δx̂ ← δx̂ + Δt · F'(x)ᵀ δx̂，与 HandwrittenTLM 严格互为转置。
// ---------------------------------------------------------------------------
class HandwrittenADM : public AdjointModel {
public:
    HandwrittenADM(dynamics::DynamicsCore& core,
                   const grid::VerticalCoordinate& vert,
                   StateReal dt);

    void propagate(const dynamics::State& x_ref, dynamics::State& dx_hat,
                   StateReal t0, StateReal t1) override;
    const char* name() const override { return "HandwrittenADM"; }

private:
    dynamics::DynamicsCore& core_;
    const grid::VerticalCoordinate& vert_;
    StateReal dt_;
};

// ---------------------------------------------------------------------------
// 内积与线性代数工具（状态空间欧氏内积）
// ---------------------------------------------------------------------------
// 状态空间的欧氏内积：⟨a, b⟩ = Σ_field Σ_i a_i b_i（覆盖 8 个字段）
StateReal state_inner_product(const dynamics::State& a,
                              const dynamics::State& b);
// axpy：y ← y + α x（状态空间）
void state_axpy(StateReal alpha, const dynamics::State& x,
                dynamics::State& y);
// 状态空间随机扰动（用于点积检验，可复现种子）
void state_random(dynamics::State& s, std::uint64_t seed, StateReal scale);

}  // namespace cubed_sph::da
