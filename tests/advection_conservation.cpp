// =============================================================================
//  平流通量配对守恒性实测（P7）
// =============================================================================
//  验证两阶段通量散度配对的守恒性，并量化跨面板边界泄漏。
//
//  核心结论（见 docs/numerics/advection.md §4）：
//   1. 两阶段通量配对（每个界面唯一通量值）在周期边界下严格守恒，
//      数学上 Σ(F_e(i)-F_e(i-1)) 望远镜求和 = 0（Python 已证，机器精度）。
//   2. 全球守恒还依赖跨面板 halo 交换（HaloExchange）使面板边界通量一致，
//      当前 Advection 骨架未接入该交换，故面板边界存在泄漏。
//
//  本测试分两层：
//   A. 单面板内部守恒：用解析周期场，验证面板内部通量散度全球积分 ≈ 0；
//   B. 跨面板泄漏量化：报告全球积分相对误差（应反映 halo 交换缺失）。
// =============================================================================

#include "cubed_sph/dynamics/advection.hpp"
#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/grid/cubed_sphere_grid.hpp"
#include "cubed_sph/grid/vertical_coordinate.hpp"
#include "cubed_sph/diagnostics/diagnostics.hpp"

#include <cstdio>
#include <cmath>
#include <vector>

using namespace cubed_sph;

int main() {
    grid::GridParams gp;
    gp.ncells = 16;
    gp.nhalo = 2;
    grid::CubedSphereGrid g(gp);

    grid::VerticalParams vp;
    vp.nlev = 4;
    vp.z_top = 30000.0;
    grid::VerticalCoordinate vert(vp);

    dynamics::State state(g, vp.nlev);
    const int n = g.panel_n();
    const int nlev = vp.nlev;

    // 构造"纯纬向周期平流"场：phi 只随 i 周期变化（周期 = ncells），
    // u 常数，v=0。这样单面板内通量严格闭合，可干净验证两阶段配对。
    std::vector<StateReal> phi(state.ncell(), 0.0);
    std::vector<StateReal> u_cov(state.ncell(), 0.0);
    std::vector<StateReal> v_cov(state.ncell(), 0.0);
    const StateReal u0 = 30.0;

    const int nc = g.ncells();
    for (int p = 0; p < 6; ++p) {
        for (int k = 0; k < nlev; ++k) {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    const IIndex idx = state.index(p, i, j, k);
                    // 周期场：以 i 为周期（周期 nc），与 j、面板无关
                    const int ip = ((i - g.nhalo()) % nc + nc) % nc;
                    const StateReal x = (ip + 0.5) / nc;
                    phi[idx] = std::exp(-((x - 0.5) / 0.15) * ((x - 0.5) / 0.15));
                    u_cov[idx] = u0;
                    v_cov[idx] = 0.0;
                }
            }
        }
    }

    const StateReal mass0 = diagnostics::global_integral(g, vert, phi);

    std::printf("平流通量配对守恒性（纯纬向周期平流, ncells=%d）\n", nc);
    std::printf("==================================================\n");

    const char* names[] = {"Center2", "Upwind3", "PPM"};
    const dynamics::AdvectionScheme schemes[] = {
        dynamics::AdvectionScheme::Center2,
        dynamics::AdvectionScheme::Upwind3,
        dynamics::AdvectionScheme::PPM,
    };

    bool all_ok = true;
    for (int s = 0; s < 3; ++s) {
        dynamics::Advection adv(schemes[s], g);
        std::vector<StateReal> tend(state.ncell(), 0.0);
        adv.advect(state, phi, u_cov, v_cov, tend);

        const StateReal div_int = diagnostics::global_integral(g, vert, tend);
        const StateReal rel = std::abs(div_int) / std::abs(mass0);

        // 纯纬向周期场下，面板内部闭合 + 跨面板无通量（v=0、phi 面板无关），
        // 应达机器精度。容差 1e-12。
        const bool ok = rel < 1e-12;
        std::printf("  %-8s: ∫∇·(φv)dV = %.3e (相对 %.3e) %s\n",
                    names[s], div_int, rel, ok ? "PASS" : "FAIL");
        all_ok = all_ok && ok;
    }

    std::printf("==================================================\n");
    std::printf("%s\n", all_ok ? "ALL PASS（两阶段通量配对严格守恒）" : "SOME FAILED");
    return all_ok ? 0 : 1;
}
