// =============================================================================
//  立方球网格度量单元测试
// =============================================================================
//  验证 gnomonic equiangular 投影的度量张量与解析值一致（提示词 P1 验收
//  标准："数值 Jacobian 与解析 Jacobian 逐点对比"）。
//
//  解析关系（Putman & Lin 2007, Eq. A1-A6）：
//    g11 = (1 + η²)/r⁴, g22 = (1 + ξ²)/r⁴, g12 = -ξη/r⁴
//    √G  = 1/r³, 其中 r² = 1 + ξ² + η²
// =============================================================================

#include <catch2/catch_test_macros.hpp>

#include "cubed_sph/grid/metric.hpp"
#include "cubed_sph/common/types.hpp"

#include <cmath>

using namespace cubed_sph;
using namespace cubed_sph::grid;

TEST_CASE("gnomonic 度量张量与解析值一致", "[grid][metric]") {
    // 选取若干非平凡 (ξ, η) 点
    const StateReal pts[][2] = {
        {0.0, 0.0}, {0.3, 0.2}, {-0.5, 0.4}, {0.6, -0.6}, {0.785398, 0.785398}};

    for (const auto& p : pts) {
        const StateReal xi = p[0];
        const StateReal eta = p[1];
        const MetricPoint m = compute_metric(xi, eta, 0);

        const StateReal r2 = 1.0 + xi * xi + eta * eta;
        const StateReal r4 = r2 * r2;

        const StateReal g11_ref = (1.0 + eta * eta) / r4;
        const StateReal g22_ref = (1.0 + xi * xi) / r4;
        const StateReal g12_ref = -xi * eta / r4;
        const StateReal sqrtG_ref = 1.0 / (r2 * std::sqrt(r2));

        REQUIRE(m.g11 == Catch::Approx(g11_ref).epsilon(1e-12));
        REQUIRE(m.g22 == Catch::Approx(g22_ref).epsilon(1e-12));
        REQUIRE(m.g12 == Catch::Approx(g12_ref).epsilon(1e-12));
        REQUIRE(m.sqrtG == Catch::Approx(sqrtG_ref).epsilon(1e-12));
    }
}

TEST_CASE("度量行列式与逆度量自洽", "[grid][metric]") {
    // det(g) 与 g·g^{-1} = I 的数值验证
    const StateReal xi = 0.4, eta = -0.3;
    const MetricPoint m = compute_metric(xi, eta, 0);

    const StateReal det = m.g11 * m.g22 - m.g12 * m.g12;

    // √G = sqrt(det)
    REQUIRE(m.sqrtG == Catch::Approx(std::sqrt(det)).epsilon(1e-12));

    // g · g^{-1} = I
    const StateReal id11 = m.g11 * m.inv_g11 + m.g12 * m.inv_g12;
    const StateReal id12 = m.g11 * m.inv_g12 + m.g12 * m.inv_g22;
    const StateReal id22 = m.g12 * m.inv_g12 + m.g22 * m.inv_g22;
    REQUIRE(id11 == Catch::Approx(1.0).epsilon(1e-12));
    REQUIRE(id12 == Catch::Approx(0.0).epsilon(1e-12));
    REQUIRE(id22 == Catch::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("协变-反变变换往返一致", "[grid][covariant]") {
    const MetricPoint m = compute_metric(0.3, 0.25, 0);
    const StateReal u_cov = 1.0, v_cov = -2.0;

    StateReal u_contra, v_contra;
    covariant_to_contravariant(m, u_cov, v_cov, u_contra, v_contra);

    StateReal u_cov2, v_cov2;
    contravariant_to_covariant(m, u_contra, v_contra, u_cov2, v_cov2);

    REQUIRE(u_cov2 == Catch::Approx(u_cov).epsilon(1e-12));
    REQUIRE(v_cov2 == Catch::Approx(v_cov).epsilon(1e-12));
}

TEST_CASE("全球面积守恒：6 面板 ∫√G dξdη = 4π", "[grid][metric]") {
    // 用梯形积分验证 6 个面板覆盖整个球面（提示词 P1 验收标准之一）。
    // 等角坐标 ξ, η ∈ [-1, 1]，单面板面积应 = 2π/3。
    const int n = 200;
    const StateReal d = 2.0 / n;
    StateReal single_panel_area = 0.0;
    for (int j = 0; j < n; ++j) {
        const StateReal eta = -1.0 + (j + 0.5) * d;
        for (int i = 0; i < n; ++i) {
            const StateReal xi = -1.0 + (i + 0.5) * d;
            single_panel_area += compute_metric(xi, eta, 0).sqrtG * d * d;
        }
    }
    const StateReal total = 6.0 * single_panel_area;
    const StateReal sphere = 4.0 * kPi;
    REQUIRE(total == Catch::Approx(sphere).epsilon(1e-3));
    REQUIRE(single_panel_area == Catch::Approx(2.0 * kPi / 3.0).epsilon(1e-3));
}
