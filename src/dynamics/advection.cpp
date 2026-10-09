// =============================================================================
//  平流实现：二阶中心 + 三阶上游 + PPM（含单调限制器）
// =============================================================================

#include "cubed_sph/dynamics/advection.hpp"

#include <algorithm>

namespace cubed_sph::dynamics {

using grid::kNumPanels;
using grid::MetricPoint;

Advection::Advection(AdvectionScheme scheme, const grid::CubedSphereGrid& grid)
    : scheme_(scheme), grid_(grid) {}

const char* Advection::scheme_name() const {
    switch (scheme_) {
        case AdvectionScheme::Center2: return "Center2";
        case AdvectionScheme::Upwind3: return "Upwind3";
        case AdvectionScheme::PPM:     return "PPM";
    }
    return "Unknown";
}

// ---------------------------------------------------------------------------
// PPM 面重构（Colella & Woodward 1984）
// ---------------------------------------------------------------------------
// 输入：一维 5 点模板 φ_{i-2}, φ_{i-1}, φ_i, φ_{i+1}, φ_{i+2}
// 返回：右界面 i+1/2 处的重构值 φ_{i+1/2}^R（迎风侧）。
//
// PPM 步骤（Colella & Woodward 1984, §1）：
//   1. 用四点插值构造单元 i 的界面值 φ_{i+1/2}：
//        φ_{i+1/2} = (7/12)(φ_i + φ_{i+1}) - (1/12)(φ_{i-1} + φ_{i+2})
//   2. 单调性限制（monotonized）：限制界面值落在相邻单元值之间，
//      抑制过冲/欠冲（保证正性）。
//
// 参考文献：
//   Colella & Woodward (1984), J. Comput. Phys. 54, 174–201
//   Lin & Rood (1996), Mon. Wea. Rev. 124, 2046–2070
// ---------------------------------------------------------------------------
inline StateReal ppm_face_value(StateReal fm2, StateReal fm1,
                                StateReal f0, StateReal fp1, StateReal fp2) {
    (void)fm2;  // 四点插值的远端点，默认形式不使用（保留接口对称性）
    // 界面 i+1/2 的四点插值（默认形式，介于 f0 与 fp1 之间）
    StateReal face = (7.0 / 12.0) * (f0 + fp1) - (1.0 / 12.0) * (fm1 + fp2);

    // 单调性限制：face 必须落在 [min(f0,fp1), max(f0,fp1)] 内
    const StateReal lo = std::min(f0, fp1);
    const StateReal hi = std::max(f0, fp1);
    face = std::max(lo, std::min(hi, face));
    return face;
}

void Advection::advect(const State& state,
                       const std::vector<StateReal>& phi,
                       const std::vector<StateReal>& u_cov,
                       const std::vector<StateReal>& v_cov,
                       std::vector<StateReal>& tendency) {
    const int n = grid_.panel_n();
    const int nlev = state.nlev();
    const StateReal dxi = 2.0 / static_cast<StateReal>(grid_.ncells());  // 等角坐标范围 [-1,1]

    // 守恒通量散度 -∇·(ρ φ v)，两阶段实现保证严格守恒：
    //   阶段一：对每个单元计算"东界面通量 F_e(i+1/2)"与"北界面通量
    //           F_n(j+1/2)"，缓存到面通量数组（每个界面一个唯一值）；
    //   阶段二：散度 = (F_e(i) - F_e(i-1) + F_n(j) - F_n(j-1)) / (√G dξ)，
    //           其中西面通量 = 西邻单元(i-1)的东面通量 F_e(i-1)，
    //           南面通量 = 南邻单元(j-1)的北面通量 F_n(j-1)。
    // 这样每个界面只被计算一次，离散严格守恒（Arakawa–Lamb 1981）。
    //
    // 面通量缓存（与 State 同布局，仅存东/北面）
    const IIndex total = state.ncell();
    std::vector<StateReal> flux_east(total, 0.0);
    std::vector<StateReal> flux_north(total, 0.0);

    // 阶段一：计算界面通量。
    // 循环范围扩到 [nhalo-1, n-nhalo]，使阶段二的西/南面通量
    // flux_east[i-1]、flux_north[j-1]（i-1 最小为 nhalo-1）有定义。
    // 面板边界处的通量闭合依赖 halo 交换（ParallelContext::HaloExchange）
    // 已填充 halo 区的 φ 与速度；本算子假设输入已含正确 halo 数据。
    for (int p = 0; p < kNumPanels; ++p) {
        for (int k = 0; k < nlev; ++k) {
            for (int j = grid_.nhalo() - 1; j < n - grid_.nhalo(); ++j) {
                for (int i = grid_.nhalo() - 1; i < n - grid_.nhalo(); ++i) {
                    const IIndex idx = state.index(p, i, j, k);

                    StateReal phi_east = 0.0, phi_north = 0.0;
                    switch (scheme_) {
                        case AdvectionScheme::Center2: {
                            phi_east = 0.5 * (phi[state.index(p, i, j, k)] +
                                              phi[state.index(p, i + 1, j, k)]);
                            phi_north = 0.5 * (phi[state.index(p, i, j, k)] +
                                               phi[state.index(p, i, j + 1, k)]);
                            break;
                        }
                        case AdvectionScheme::Upwind3: {
                            const StateReal u_e = u_cov[idx];
                            const StateReal v_n = v_cov[idx];
                            if (u_e >= 0.0) {
                                phi_east = (2.0 * phi[state.index(p, i + 1, j, k)] +
                                            5.0 * phi[state.index(p, i, j, k)] -
                                            phi[state.index(p, i - 1, j, k)]) / 6.0;
                            } else {
                                phi_east = (2.0 * phi[state.index(p, i, j, k)] +
                                            5.0 * phi[state.index(p, i + 1, j, k)] -
                                            phi[state.index(p, i + 2, j, k)]) / 6.0;
                            }
                            if (v_n >= 0.0) {
                                phi_north = (2.0 * phi[state.index(p, i, j + 1, k)] +
                                             5.0 * phi[state.index(p, i, j, k)] -
                                             phi[state.index(p, i, j - 1, k)]) / 6.0;
                            } else {
                                phi_north = (2.0 * phi[state.index(p, i, j, k)] +
                                             5.0 * phi[state.index(p, i, j + 1, k)] -
                                             phi[state.index(p, i, j + 2, k)]) / 6.0;
                            }
                            break;
                        }
                        case AdvectionScheme::PPM: {
                            phi_east = ppm_face_value(
                                phi[state.index(p, i - 1, j, k)],
                                phi[state.index(p, i, j, k)],
                                phi[state.index(p, i, j, k)],
                                phi[state.index(p, i + 1, j, k)],
                                phi[state.index(p, i + 2, j, k)]);
                            phi_north = ppm_face_value(
                                phi[state.index(p, i, j - 1, k)],
                                phi[state.index(p, i, j, k)],
                                phi[state.index(p, i, j, k)],
                                phi[state.index(p, i, j + 1, k)],
                                phi[state.index(p, i, j + 2, k)]);
                            break;
                        }
                    }
                    // 东/北界面通量（界面速度取本单元面速度 u_cov/v_cov）
                    flux_east[idx]  = phi_east * u_cov[idx];
                    flux_north[idx] = phi_north * v_cov[idx];
                }
            }
        }
    }

    // 阶段二：通量散度（严格守恒配对）
    for (int p = 0; p < kNumPanels; ++p) {
        for (int k = 0; k < nlev; ++k) {
            for (int j = grid_.nhalo(); j < n - grid_.nhalo(); ++j) {
                for (int i = grid_.nhalo(); i < n - grid_.nhalo(); ++i) {
                    const IIndex idx = state.index(p, i, j, k);
                    const MetricPoint& m = grid_.metric(p, i, j);
                    // 西面通量 = 西邻单元的东面通量；南面通量 = 南邻的北面通量
                    const StateReal fe = flux_east[idx];
                    const StateReal fw = flux_east[state.index(p, i - 1, j, k)];
                    const StateReal fn = flux_north[idx];
                    const StateReal fs = flux_north[state.index(p, i, j - 1, k)];
                    tendency[idx] += -(fe - fw + fn - fs) / (m.sqrtG * dxi);
                }
            }
        }
    }
}

}  // namespace cubed_sph::dynamics
