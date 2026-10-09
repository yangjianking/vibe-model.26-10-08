// =============================================================================
//  P7 自检：守恒诊断 + 标准算例逻辑验证
// =============================================================================
//  验证：
//   1. 全球积分与 4π 面积守恒（度量加权正确）
//   2. 等温静力平衡态的总能量/质量诊断合理性
//   3. Williamson TC2 解析解（钟型场归一化 + 12 天周期平移）
// =============================================================================

#include "cubed_sph/diagnostics/diagnostics.hpp"
#include "cubed_sph/testcases/test_cases.hpp"
#include "cubed_sph/dynamics/thermo.hpp"
#include "cubed_sph/grid/vertical_coordinate.hpp"

#include <cstdio>
#include <cmath>
#include <vector>

using namespace cubed_sph;

int main() {
    // 网格（ncells=16, nhalo=2, nlev=4）
    grid::GridParams gp;
    gp.ncells = 16;
    gp.nhalo = 2;
    grid::CubedSphereGrid g(gp);

    grid::VerticalParams vp;
    vp.nlev = 16;
    vp.z_top = 30000.0;
    grid::VerticalCoordinate vert(vp);

    dynamics::State state(g, vp.nlev);

    // ---- 测试 1：全球积分 = 单位场积分 ≈ 4π R² × 层厚总和 ----
    // 用常密度 1 场，全球积分应 ≈ 4π R² · Σ dz（地球表面积 × 总高度）
    std::vector<StateReal> ones(state.ncell(), 1.0);
    StateReal vol = diagnostics::global_integral(g, vert, ones);
    const StateReal R = g.radius();
    const StateReal area = 4.0 * 3.141592653589793 * R * R;
    const StateReal total_dz = vert.z_top();  // z_half 从 0 到 z_top
    const StateReal expected = area * total_dz;
    StateReal rel1 = std::abs(vol - expected) / expected;
    std::printf("[1] 全球积分(单位场): vol=%.6e expected=%.6e rel=%.3e %s\n",
                vol, expected, rel1, rel1 < 2e-3 ? "PASS" : "FAIL");

    // ---- 测试 2：等温静力平衡态质量诊断 ----
    // 用 TC5 初始化（等温静力平衡），计算总质量应接近解析值
    testcases::WilliamsonTC5 tc5;
    tc5.initialize(state, vert);
    diagnostics::ConservationDiagnostics diag(g, vert);
    auto snap = diag.compute(state);
    // 解析总质量 = ∫ ρ0 exp(-z/H) dV = ρ0 H (1 - exp(-z_top/H)) * area
    const StateReal T0 = 288.0;
    const StateReal H = 287.05 * T0 / 9.80665;
    const StateReal rho0 = 100000.0 / (287.05 * T0);
    const StateReal expected_mass = rho0 * H * (1.0 - std::exp(-vert.z_top() / H)) * area;
    StateReal rel2 = std::abs(snap.total_mass - expected_mass) / expected_mass;
    std::printf("[2] 静力平衡态总质量: mass=%.6e expected=%.6e rel=%.3e %s\n",
                snap.total_mass, expected_mass, rel2, rel2 < 2e-2 ? "PASS" : "FAIL");

    // ---- 测试 3：Williamson TC2 解析解归一化 ----
    testcases::WilliamsonTC2 tc2;
    // 钟型场峰值应为 1.0（在中心），远端为 0
    StateReal peak = tc2.analytic(tc2.lambda_c, tc2.phi_c, 0.0);
    StateReal far = tc2.analytic(tc2.lambda_c + 3.141592653589793,
                                 tc2.phi_c, 0.0);  // 对跖点
    std::printf("[3] TC2 解析解: peak=%.6f far=%.6f %s\n",
                peak, far,
                (std::abs(peak - 1.0) < 1e-12 && far < 1e-12) ? "PASS" : "FAIL");

    // ---- 测试 4：TC2 12 天周期平移（解析解自洽）----
    // 12 天后场回到原点：h(lon,lat,12天) == h(lon,lat,0)
    const StateReal T12 = 12.0 * 86400.0;
    StateReal h0 = tc2.analytic(0.3, 0.2, 0.0);
    StateReal h12 = tc2.analytic(0.3, 0.2, T12);
    StateReal rel4 = std::abs(h0 - h12);
    std::printf("[4] TC2 12天周期: h(0)=%.6f h(12d)=%.6f diff=%.3e %s\n",
                h0, h12, rel4, rel4 < 1e-6 ? "PASS" : "FAIL");

    // ---- 测试 5：Held-Suarez 平衡温度剖面单调性 ----
    StateReal T_eq_lo = testcases::HeldSuarez::equilibrium_temperature(0.0, 0.0);
    StateReal T_eq_hi = testcases::HeldSuarez::equilibrium_temperature(30000.0, 0.0);
    StateReal T_eq_pole = testcases::HeldSuarez::equilibrium_temperature(0.0, 1.4);
    std::printf("[5] Held-Suarez: T(赤道低层)=%.1f T(高层)=%.1f T(极地)=%.1f %s\n",
                T_eq_lo, T_eq_hi, T_eq_pole,
                (T_eq_lo > T_eq_pole) ? "PASS" : "FAIL");

    // ---- 测试 6：瑞利摩擦系数边界层性质 ----
    // σ = 1 - z/z_top，低层 σ→1（z 小）应有摩擦，高层 σ→0（z 大）无摩擦
    StateReal kv_surf = testcases::HeldSuarez::rayleigh_friction(3000.0);   // σ≈0.9
    StateReal kv_top = testcases::HeldSuarez::rayleigh_friction(24000.0);   // σ≈0.2
    std::printf("[6] Held-Suarez 摩擦: kv(σ≈0.9)=%.2e kv(σ≈0.2)=%.2e %s\n",
                kv_surf, kv_top, (kv_surf > kv_top) ? "PASS" : "FAIL");

    bool ok = (rel1 < 2e-3) && (rel2 < 2e-2) &&
              (std::abs(peak - 1.0) < 1e-12 && far < 1e-12) &&
              (rel4 < 1e-6) && (T_eq_lo > T_eq_pole) && (kv_surf > kv_top);
    std::printf("\n=== %s ===\n", ok ? "ALL PASS" : "SOME FAILED");
    return ok ? 0 : 1;
}
