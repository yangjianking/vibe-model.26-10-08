// =============================================================================
//  守恒诊断（Conservation Diagnostics）（P7 检验评估）
// =============================================================================
//  计算全球守恒量及其随时间的漂移，作为模式正确性与长期积分稳定性的
//  核心检验手段。非静力可压缩 Euler 方程（绝热、无摩擦、无地形）应精确
//  守恒以下量（离散意义下，Arakawa–Lamb 1981 配对保证）：
//
//    1. 总质量        M  = ∫ ρ dV
//    2. 总位温        Θ  = ∫ ρθ dV      （绝热无源时守恒）
//    3. 总能量        E  = ∫ (½ρ|v|² + ρ c_v T + ρ Φ) dV
//    4. 总角动量      A  = ∫ ρ (u r cosφ + Ω r² cos²φ) dV
//
//  全球积分在立方球上为 6 个面板的积分之和，每个面板用 gnomonic 度量
//  √G 加权。采用 Kahan/长双精度补偿求和，避免大网格求和误差污染守恒诊断
//  （提示词 P3 第 2 条：禁止朴素求和）。
//
//  参考文献：
//    Arakawa & Lamb (1981), Mon. Wea. Rev. 109, 18–36
//    Williamson et al. (1992), J. Comput. Phys. 102, 211–224（标准算例）
//    Thuburn (2008), J. Comput. Phys. 227, 3715–3730（离散守恒）
// =============================================================================

#pragma once

#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/grid/vertical_coordinate.hpp"
#include "cubed_sph/grid/metric.hpp"

#include <vector>

namespace cubed_sph::diagnostics {

// ---------------------------------------------------------------------------
// 全球守恒量快照
// ---------------------------------------------------------------------------
struct ConservationSnapshot {
    StateReal total_mass = 0.0;      // ∫ ρ dV
    StateReal total_theta = 0.0;     // ∫ ρθ dV
    StateReal total_energy = 0.0;    // ∫ (½ρ|v|² + ρc_vT + ρΦ) dV
    StateReal total_angular_momentum = 0.0;  // ∫ ρ(u r cosφ + Ω r² cos²φ) dV

    StateReal kinetic_energy = 0.0;  // ∫ ½ρ|v|² dV（能量分解）
    StateReal internal_energy = 0.0; // ∫ ρ c_v T dV
    StateReal potential_energy = 0.0;// ∫ ρ Φ dV
};

// ---------------------------------------------------------------------------
// 守恒诊断器
// ---------------------------------------------------------------------------
// 绑定网格与垂直坐标，提供全球积分与守恒误差（相对初值）的计算。
// ---------------------------------------------------------------------------
class ConservationDiagnostics {
public:
    ConservationDiagnostics(const grid::CubedSphereGrid& grid,
                            const grid::VerticalCoordinate& vert);

    // 计算当前状态的全球守恒量快照
    ConservationSnapshot compute(const dynamics::State& state) const;

    // 计算相对初值的守恒误差（相对误差 = |X(t) - X(0)| / |X(0)|）
    // 返回各量的相对漂移。
    static std::vector<StateReal> relative_drift(
        const ConservationSnapshot& initial,
        const ConservationSnapshot& current);

    // 网格与垂直坐标访问
    const grid::CubedSphereGrid& grid() const { return grid_; }
    const grid::VerticalCoordinate& vert() const { return vert_; }

private:
    const grid::CubedSphereGrid& grid_;
    const grid::VerticalCoordinate& vert_;
};

// ---------------------------------------------------------------------------
// 全球积分工具（供算例与诊断共用）
// ---------------------------------------------------------------------------
// 对 6 面板做 gnomonic 度量加权的全球积分 ∫ f √G dξ dη dη_vert，
// 采用长双精度补偿求和。f 为单元中心标量场（与 State 同布局）。
// ---------------------------------------------------------------------------
StateReal global_integral(const grid::CubedSphereGrid& grid,
                          const grid::VerticalCoordinate& vert,
                          const std::vector<StateReal>& field);

}  // namespace cubed_sph::diagnostics
