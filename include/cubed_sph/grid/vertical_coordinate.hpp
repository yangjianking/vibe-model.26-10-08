// =============================================================================
//  垂直坐标子系统（Vertical Coordinate）
// =============================================================================
//  垂直离散采用地形跟随混合坐标（可选纯高度坐标）。本模块提供：
//    - 垂直层级（level）定义与交错（Charney–Phillips / Lorenz）
//    - 地形跟随坐标变换 η(x, y, z) 及其度量项（∂z/∂η、∂η/∂z、Jacobian）
//    - 地形平滑（Schär et al. 2002）与陡坡限制（Janjić 1984 判据）
//
//  坐标定义（地形跟随混合 σ–p / 高度坐标，Simmons & Burridge 1981；
//  采用高度型地形跟随坐标，Schär et al. 2002 平滑地形跟随）：
//
//    z(η) = z_top * η + h_surf * b(η)
//
//  其中：
//    η      ∈ [0, 1]：归一化垂直坐标（0=地面，1=模式顶）
//    z_top  ：模式顶高度
//    h_surf ：地形高度（已平滑）
//    b(η)   ：地形衰减函数（Schär 2002 的 b(η) = sinh{(η_top - η)/s}
//             / sinh{η_top/s} 形式，或用多项式），在地面 b(0)=1，模式顶 b(1)=0。
//
//  度量项（垂直 Jacobian）：
//    J = ∂z/∂η = z_top + h_surf * db/dη
//    ∂η/∂z = 1 / J
//
//  地形斜率限制（Janjić 1984 类判据）：限制 ∂h/∂x、∂h/∂y 与 Δh/Δz 的比值，
//  避免陡峭地形导致坐标面交叉与数值不稳定。
//
//  参考文献：
//    Phillips (1956), J. Meteorol. 14, 184–185（σ 坐标）
//    Simmons & Burridge (1981), Mon. Wea. Rev. 109, 758–766（混合坐标）
//    Schär et al. (2002), Mon. Wea. Rev. 130, 2459–2480（平滑地形跟随）
//    Janjić (1984) [TO-VERIFY]（地形斜率限制）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"

#include <vector>

namespace cubed_sph::grid {

// ---------------------------------------------------------------------------
// 垂直交错方案
// ---------------------------------------------------------------------------
enum class VerticalStagger {
    Lorenz,           // 速度与温度/气压同位（经典 Lorenz 交错）
    CharneyPhillips,  // 温度与速度交错（避免虚假的静力不平衡）
};

// ---------------------------------------------------------------------------
// 垂直坐标类型
// ---------------------------------------------------------------------------
enum class VerticalCoordType {
    TerrainFollowing,  // 地形跟随混合坐标（默认）
    PureHeight,        // 纯高度坐标（z 坐标）
};

// ---------------------------------------------------------------------------
// 垂直坐标参数（由配置注入）
// ---------------------------------------------------------------------------
struct VerticalParams {
    VerticalCoordType type = VerticalCoordType::TerrainFollowing;
    VerticalStagger stagger = VerticalStagger::CharneyPhillips;
    int nlev = 60;            // 垂直层数
    StateReal z_top = 30000.0; // 模式顶高度 [m]
    StateReal s = 5.0;         // Schär 平滑参数 s（地形衰减尺度，无量纲）
    StateReal max_slope = 0.5; // 地形斜率上限（Janjić 判据，无量纲 Δh/Δx）
};

// ---------------------------------------------------------------------------
// 垂直坐标
// ---------------------------------------------------------------------------
// 预计算并缓存：层级高度 z_k、η 坐标、垂直 Jacobian J_k = ∂z/∂η、∂η/∂z、
// 地形衰减函数 b(η) 及其导数。供动力核心与 I/O 共享。
// ---------------------------------------------------------------------------
class VerticalCoordinate {
public:
    explicit VerticalCoordinate(const VerticalParams& params);

    // ---- 维度与参数 ----
    int nlev() const { return nlev_; }
    StateReal z_top() const { return z_top_; }
    VerticalCoordType type() const { return type_; }
    VerticalStagger stagger() const { return stagger_; }

    // ---- 层级高度 ----
    // 主层（full level）高度 z_k（温度/湿度所在层）
    StateReal z_full(int k) const { return z_full_[k]; }
    // 半层（half level / 界面）高度 z_{k+1/2}（垂直速度/通量所在层）
    StateReal z_half(int k) const { return z_half_[k]; }

    // ---- 度量项 ----
    // 主层垂直 Jacobian J_k = (∂z/∂η)_k
    StateReal jacobian(int k) const { return jacobian_[k]; }
    // 半层 ∂η/∂z
    StateReal d_eta_dz(int k) const { return deta_dz_[k]; }

    // ---- 地形 ----
    // 设置平滑后的地形高度（在动力学初始化时由地形数据填充）
    void set_surface_height(StateReal h) { h_surf_ = h; }
    StateReal surface_height() const { return h_surf_; }

    // 重建度量（地形高度改变后调用）
    void rebuild();

    // 地形衰减函数 b(η) 及其导数 db/dη（Schär 2002）
    static StateReal decay(StateReal eta, StateReal s);
    static StateReal decay_deriv(StateReal eta, StateReal s);

private:
    VerticalParams params_;
    int nlev_;
    StateReal z_top_;
    VerticalCoordType type_;
    VerticalStagger stagger_;
    StateReal h_surf_ = 0.0;

    std::vector<StateReal> z_full_;   // 主层高度 [m]
    std::vector<StateReal> z_half_;   // 半层高度 [m]（nlev+1 个）
    std::vector<StateReal> jacobian_; // 主层 J = ∂z/∂η
    std::vector<StateReal> deta_dz_;  // 半层 ∂η/∂z

    void build_levels();
};

}  // namespace cubed_sph::grid
