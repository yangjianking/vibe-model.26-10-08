// =============================================================================
//  地形处理实现
// =============================================================================

#include "cubed_sph/grid/topography.hpp"

#include <cmath>
#include <algorithm>

namespace cubed_sph::grid {

namespace {
// 面板内一维扁平索引
inline IIndex flat(const CubedSphereGrid& g, int p, int i, int j) {
    return static_cast<IIndex>(p) * g.panel_n() * g.panel_n() +
           static_cast<IIndex>(j) * g.panel_n() + i;
}
}  // namespace

void Topography::smooth(const CubedSphereGrid& grid,
                        std::vector<StateReal>& h,
                        int n_smooth_passes,
                        StateReal smoothing_factor) {
    const int n = grid.panel_n();
    // 多次 5 点拉普拉斯平滑：h_new = h + factor * (∇²h)/4
    // 平滑因子 ∈ [0,1]，越大平滑越强。
    std::vector<StateReal> h_new(h.size());
    for (int pass = 0; pass < n_smooth_passes; ++pass) {
        for (int p = 0; p < kNumPanels; ++p) {
            for (int j = grid.nhalo(); j < n - grid.nhalo(); ++j) {
                for (int i = grid.nhalo(); i < n - grid.nhalo(); ++i) {
                    const StateReal c = h[flat(grid, p, i, j)];
                    const StateReal e = h[flat(grid, p, i + 1, j)];
                    const StateReal w = h[flat(grid, p, i - 1, j)];
                    const StateReal n0 = h[flat(grid, p, i, j + 1)];
                    const StateReal s0 = h[flat(grid, p, i, j - 1)];
                    const StateReal lap = (e + w + n0 + s0) / 4.0 - c;
                    h_new[flat(grid, p, i, j)] = c + smoothing_factor * lap;
                }
            }
        }
        std::swap(h, h_new);
    }
}

void Topography::limit_slope(const CubedSphereGrid& grid,
                             std::vector<StateReal>& h,
                             StateReal max_slope,
                             StateReal dx_meters) {
    const int n = grid.panel_n();
    // Janjić 判据：|Δh| ≤ max_slope * Δx，超出部分向低处削峰（限制坡度）。
    for (int p = 0; p < kNumPanels; ++p) {
        for (int j = grid.nhalo(); j < n - grid.nhalo(); ++j) {
            for (int i = grid.nhalo(); i < n - grid.nhalo(); ++i) {
                const StateReal hmax = max_slope * dx_meters;
                const StateReal c = h[flat(grid, p, i, j)];
                const StateReal e = h[flat(grid, p, i + 1, j)];
                const StateReal n0 = h[flat(grid, p, i, j + 1)];
                // 限制东/北向相邻点的高差
                if (e - c > hmax) h[flat(grid, p, i + 1, j)] = c + hmax;
                if (c - e > hmax) h[flat(grid, p, i + 1, j)] = c - hmax;
                if (n0 - c > hmax) h[flat(grid, p, i, j + 1)] = c + hmax;
                if (c - n0 > hmax) h[flat(grid, p, i, j + 1)] = c - hmax;
            }
        }
    }
}

void Topography::gaussian_mountain(const CubedSphereGrid& grid,
                                   std::vector<StateReal>& h,
                                   StateReal lon0, StateReal lat0,
                                   StateReal half_width, StateReal height) {
    const int n = grid.panel_n();
    const StateReal r = grid.radius();
    for (int p = 0; p < kNumPanels; ++p) {
        for (int j = grid.nhalo(); j < n - grid.nhalo(); ++j) {
            for (int i = grid.nhalo(); i < n - grid.nhalo(); ++i) {
                const MetricPoint& m = grid.metric(p, i, j);
                // 球面大圆距离近似（用经纬度差，小尺度地形下足够）
                const StateReal dlon = m.lon - lon0;
                const StateReal dlat = m.lat - lat0;
                const StateReal coslat = std::cos(lat0);
                const StateReal dist = r * std::sqrt(dlon * dlon * coslat * coslat +
                                                     dlat * dlat);
                const StateReal r2 = (dist / half_width) * (dist / half_width);
                h[flat(grid, p, i, j)] = height * std::exp(-r2);
            }
        }
    }
}

}  // namespace cubed_sph::grid
