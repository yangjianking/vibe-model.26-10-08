// =============================================================================
//  面板旋转矩阵构造（离线预计算）
// =============================================================================
//  基于面板局部基 (e1, e2) 的朝向，构造面板间切向旋转矩阵。
//  旋转矩阵元素为两面板基向量的内积（见 panel_rotation.hpp 的推导）。
// =============================================================================

#include "cubed_sph/grid/panel_rotation.hpp"

#include <cmath>
#include <array>

namespace cubed_sph::grid {

namespace {

// 面板局部基向量（与 cubed_sphere_grid.cpp 中 kPanelE1/kPanelE2 一致）
constexpr std::array<std::array<StateReal, 3>, 6> kE1 = {{
    {{ 0.0, 1.0, 0.0}}, {{-1.0, 0.0, 0.0}}, {{ 0.0, -1.0, 0.0}},
    {{ 1.0, 0.0, 0.0}}, {{ 1.0, 0.0, 0.0}},  {{ 1.0, 0.0, 0.0}},
}};

constexpr std::array<std::array<StateReal, 3>, 6> kE2 = {{
    {{ 0.0, 0.0, 1.0}}, {{ 0.0, 0.0, 1.0}}, {{ 0.0, 0.0, 1.0}},
    {{ 0.0, 0.0, 1.0}}, {{ 0.0, 1.0, 0.0}}, {{ 0.0, -1.0, 0.0}},
}};

inline StateReal dot3(const std::array<StateReal, 3>& a,
                      const std::array<StateReal, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

}  // namespace

PanelRotation panel_rotation(PanelId src, PanelId dst) {
    const int s = static_cast<int>(src);
    const int d = static_cast<int>(dst);
    PanelRotation R;
    // R = [[e1_s·e1_d, e2_s·e1_d],
    //      [e1_s·e2_d, e2_s·e2_d]]
    R.r11 = dot3(kE1[s], kE1[d]);
    R.r12 = dot3(kE2[s], kE1[d]);
    R.r21 = dot3(kE1[s], kE2[d]);
    R.r22 = dot3(kE2[s], kE2[d]);
    return R;
}

}  // namespace cubed_sph::grid
