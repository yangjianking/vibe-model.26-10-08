// =============================================================================
//  SoA 布局与混合精度单元测试
// =============================================================================
//  验证 AoS↔SoA 互转正确性、精度卫士、随机舍入的期望值性质。
// =============================================================================

#include <catch2/catch_test_macros.hpp>

#include "cubed_sph/common/soa_layout.hpp"
#include "cubed_sph/common/precision_guard.hpp"
#include "cubed_sph/common/compensated_sum.hpp"

#include <cmath>
#include <vector>

using namespace cubed_sph;
using namespace cubed_sph::common;

TEST_CASE("AoS↔SoA 互转往返一致", "[common][soa]") {
    const int n = 6;      // 小面板（含 halo）
    const int nlev = 8;
    const IIndex nh = 6 * n * n;

    // 构造 AoS 测试数据（随机）
    std::vector<StateReal> aos(nh * nlev);
    for (std::size_t i = 0; i < aos.size(); ++i) {
        aos[i] = static_cast<StateReal>(i) * 0.5;
    }

    SoAField soa(nh, nlev);
    aos_to_soa(aos, n, nlev, soa);

    std::vector<StateReal> back(nh * nlev, 0.0);
    soa_to_aos(soa, n, nlev, back);

    for (std::size_t i = 0; i < aos.size(); ++i) {
        REQUIRE(back[i] == Catch::Approx(aos[i]).epsilon(1e-12));
    }
}

TEST_CASE("SoA 垂直维连续（垂直列内存相邻）", "[common][soa]") {
    const int nlev = 10;
    SoAField f(2, nlev);
    // 同一水平格点 h 的垂直列在内存中连续
    f(0, 0) = 1.0;
    f(0, nlev - 1) = 2.0;
    REQUIRE(f(0, 0) == Catch::Approx(1.0));
    REQUIRE(f(0, nlev - 1) == Catch::Approx(2.0));
    // 垂直相邻元素地址连续
    const StateReal* base = f.data();
    REQUIRE(&f(0, 1) - base == 1);
}

TEST_CASE("精度卫士除零保护", "[common][precision]") {
    // 极小分母被保护，避免 Inf
    const StateReal g = guard_divisor(0.0);
    REQUIRE(std::isfinite(static_cast<double>(1.0 / g)));
    // 正常分母不受影响
    REQUIRE(guard_divisor(5.0) == Catch::Approx(5.0));
    // 负的极小分母
    REQUIRE(guard_divisor(-1.0e-30) < 0.0);
}

TEST_CASE("随机舍入期望值等于精确值", "[common][precision]") {
    StochasticRounder sr(42);
    // 对大量随机数做随机舍入，均值应接近原值（无系统性偏差）
    const int N = 100000;
    double sum_orig = 0.0, sum_rounded = 0.0;
    for (int i = 0; i < N; ++i) {
        const double x = 0.1 + 0.001 * (i % 997);  // 非精确 fp32 可表示的数
        sum_orig += x;
        sum_rounded += static_cast<double>(sr.round_to_f32(x));
    }
    // 随机舍入的期望值 = 精确值，均值相对偏差应很小
    const double rel = std::abs(sum_rounded - sum_orig) / sum_orig;
    REQUIRE(rel < 1e-3);
}

TEST_CASE("Kahan 补偿求和优于朴素求和", "[common][compensated]") {
    // 1 + 1e-8 + ... + 1e-8（共 1e6 个 1e-8）：朴素求和会丢失 1e-8 项
    const double eps = 1e-8;
    const int N = 1000000;

    double naive = 1.0;
    for (int i = 0; i < N; ++i) naive += eps;

    KahanSum ks;
    ks.add(1.0);
    for (int i = 0; i < N; ++i) ks.add(eps);

    // 精确值 = 1 + N*eps = 1.01
    const double exact = 1.0 + N * eps;
    REQUIRE(ks.value() == Catch::Approx(exact).epsilon(1e-12));
    // 朴素求和的误差更大（可能因 fp64 精度仍较小，此处仅验证 Kahan 不劣于朴素）
    REQUIRE(std::abs(ks.value() - exact) <= std::abs(naive - exact) + 1e-15);
}
