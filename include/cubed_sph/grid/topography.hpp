// =============================================================================
//  地形处理（Topography）：平滑与陡坡限制
// =============================================================================
//  地形跟随坐标要求地形足够平滑，避免坐标面交叉与虚假地形强迫。本模块提供：
//    - 地形平滑（高斯/拉普拉斯平滑，Schär et al. 2002 的 scale-selective 平滑）
//    - 陡坡限制（Janjić 1984 类判据：限制 Δh/Δx，超出部分削峰）
//
//  平滑动机：地形跟随坐标在高分辨率下的数值稳定性要求地形谱的短波分量被
//  有效抑制（Schär et al. 2002）。限制器保证地形斜率不超过配置上限。
//
//  参考文献：
//    Schär et al. (2002), Mon. Wea. Rev. 130, 2459–2480
//    Janjić (1984) [TO-VERIFY]（地形斜率限制）
// =============================================================================

#pragma once

#include "cubed_sph/grid/cubed_sphere_grid.hpp"
#include "cubed_sph/common/types.hpp"

#include <vector>

namespace cubed_sph::grid {

// ---------------------------------------------------------------------------
// 地形处理
// ---------------------------------------------------------------------------
class Topography {
public:
    // 对面板上的原始地形高度做平滑（多次拉普拉斯平滑）
    // 输入/输出：与网格同布局（panel × n × n）的地形高度 [m]
    static void smooth(const CubedSphereGrid& grid,
                       std::vector<StateReal>& h,
                       int n_smooth_passes,
                       StateReal smoothing_factor);

    // 陡坡限制（Janjić 判据）：对每个面板限制 |Δh/Δx| ≤ max_slope
    static void limit_slope(const CubedSphereGrid& grid,
                            std::vector<StateReal>& h,
                            StateReal max_slope,
                            StateReal dx_meters);

    // 理想化地形（供后续阶段标准算例使用）
    //  - gaussian_mountain：孤立高斯山（Williamson Test 5 风格）
    static void gaussian_mountain(const CubedSphereGrid& grid,
                                  std::vector<StateReal>& h,
                                  StateReal lon0, StateReal lat0,
                                  StateReal half_width, StateReal height);
};

}  // namespace cubed_sph::grid
