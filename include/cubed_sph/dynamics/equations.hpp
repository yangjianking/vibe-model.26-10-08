// =============================================================================
//  动力方程（Equations）—— 非静力可压缩 Euler 方程的通量形式
// =============================================================================
//  本模块提供动力核心的右端项（倾向）计算，即"连续形式 → 离散形式"的
//  落地。完整推导见 docs/theory/dynamics_equations.md，此处给出实现所依据
//  的离散公式与守恒配对。
//
//  预后变量：ρ, ρu, ρv, ρw, ρθ（见 state.hpp）。
//
//  连续方程（质量守恒，通量形式）：
//    ∂ρ/∂t = -∇·(ρ v)
//
//  动量方程（守恒形式，含科氏力与气压梯度）：
//    ∂(ρu)/∂t = -∇·(ρu v) - ∂p/∂ξ + ρ F_coriolis + ...   （水平，协变）
//    ∂(ρw)/∂t = -∇·(ρw v) - ∂p/∂z - ρ g                     （垂直）
//
//  热力学方程（位温密度守恒）：
//    ∂(ρθ)/∂t = -∇·(ρθ v) + 源项（物理过程）
//
//  守恒配对（Arakawa–Lamb 1981）：散度算子与压力梯度算子采用互为伴随的
//  离散形式，保证离散总能量守恒。这里将气压梯度写成 Exner 压力梯度与
//  位温的耦合形式（见 thermo.hpp 与 docs/numerics/pressure_gradient.md）。
//
//  参考文献：
//    Arakawa & Lamb (1981), Mon. Wea. Rev. 109, 18–36（能量守恒配对）
//    Lin & Rood (1996), Mon. Wea. Rev. 124, 2046–2070（通量形式平流）
//    Harris et al. (2021), GFDL TM GFDL2021001（FV³ 方程集）
//    Wood et al. (2014), QJRMS 140（ENDGame 半隐式方程）
// =============================================================================

#pragma once

#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/dynamics/thermo.hpp"
#include "cubed_sph/grid/differential_operators.hpp"
#include "cubed_sph/grid/vertical_coordinate.hpp"
#include "cubed_sph/common/config.hpp"

#include <vector>

namespace cubed_sph::dynamics {

// ---------------------------------------------------------------------------
// 倾向容器（与 State 同布局）
// ---------------------------------------------------------------------------
struct Tendency {
    explicit Tendency(const grid::CubedSphereGrid& g, int nlev);

    std::vector<StateReal> drho;        // ∂ρ/∂t
    std::vector<StateReal> drho_u;      // ∂(ρu)/∂t
    std::vector<StateReal> drho_v;      // ∂(ρv)/∂t
    std::vector<StateReal> drho_w;      // ∂(ρw)/∂t
    std::vector<StateReal> drho_theta;  // ∂(ρθ)/∂t

    // 清零（每个时间步开始时调用）
    void zero();
};

// ---------------------------------------------------------------------------
// 动力核心：计算右端项（显式部分）
// ---------------------------------------------------------------------------
// 职责：
//   1. 质量/热力学平流（通量形式，Lin–Rood 1996 风格，占位接口）
//   2. 气压梯度力（Exner 压力梯度 + 位温，Arakawa–Lamb 配对）
//   3. 科氏力与曲率项（协变形式，见 coriolis.hpp）
//   4. 垂直方向（重力、垂直平流、垂直气压梯度）
//
// 半隐式处理：快波线性项（气压梯度与散度耦合、浮力与垂直速度耦合）在
// timeint/ 中用 Crank–Nicolson 隐式平均，本函数只负责显式残差部分。
// ---------------------------------------------------------------------------
class DynamicsCore {
public:
    DynamicsCore(const grid::CubedSphereGrid& grid,
                 const grid::VerticalCoordinate& vert,
                 int nlev,
                 const Config& cfg);

    // 计算显式右端项（不含半隐式隐式项）
    void compute_explicit(const State& state, Tendency& tendency);

    // 应用科氏力与曲率项（协变动量方程的源项）
    void apply_coriolis(const State& state, Tendency& tendency);

    // 应用垂直方向项（重力、垂直气压梯度、垂直平流）
    void apply_vertical(const State& state, Tendency& tendency);

private:
    const grid::CubedSphereGrid& grid_;
    const grid::VerticalCoordinate& vert_;
    int nlev_;
    StateReal coriolis_f_;   // 科氏参数参考值（2Ω sin φ₀，可配）
};

}  // namespace cubed_sph::dynamics
