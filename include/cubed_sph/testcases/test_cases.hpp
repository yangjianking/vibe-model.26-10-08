// =============================================================================
//  标准算例（Test Cases）（P7 检验评估）
// =============================================================================
//  Williamson et al. (1992) 提出的浅水标准算例集（TC1–TC7）是球面全球
//  大气模式正确性验证的黄金标准。此外，Jablonowski & Williamson (2006)
//  的斜压波算例与 Held & Suarez (1994) 的气候态强迫用于三维模式。
//
//  本模块提供：
//    - Williamson TC2：稳态定常平流（有解析解，可做误差收敛阶验证）
//    - Williamson TC5：地形上的罗斯贝波（zonal flow over mountain）
//    - Jablonowski–Williamson：斜压波（三维，含垂直结构）
//    - Held–Suarez：气候态温度强迫/瑞利摩擦（三维长期积分）
//
//  各算例提供初始条件（填充 State）与解析解（TC2 有，供误差计算）。
//
//  参考文献：
//    Williamson et al. (1992), J. Comput. Phys. 102, 211–224
//    Jablonowski & Williamson (2006), QJRMS 132, 2943–2975
//    Held & Suarez (1994), Bull. Amer. Meteor. Soc. 75, 1825–1830
// =============================================================================

#pragma once

#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/grid/cubed_sphere_grid.hpp"
#include "cubed_sph/grid/vertical_coordinate.hpp"

namespace cubed_sph::testcases {

// ---------------------------------------------------------------------------
// 算例类型枚举
// ---------------------------------------------------------------------------
enum class TestCaseType {
    WilliamsonTC2,          // 稳态定常平流（解析解）
    WilliamsonTC5,          // 地形罗斯贝波
    JablonowskiWilliamson,  // 斜压波（三维）
    HeldSuarez,             // 气候态强迫
};

// ---------------------------------------------------------------------------
// 算例基类：初始化状态（初始条件）
// ---------------------------------------------------------------------------
class TestCase {
public:
    virtual ~TestCase() = default;

    // 填充初始条件到 state（ρ、ρθ、动量密度等）
    virtual void initialize(dynamics::State& state,
                            const grid::VerticalCoordinate& vert) = 0;

    virtual const char* name() const = 0;
};

// ---------------------------------------------------------------------------
// Williamson TC2：稳态定常平流（有解析解）
// ---------------------------------------------------------------------------
// 一个余弦钟型标量场被恒定角速度的定常风场平流，无变形。解析解为钟型场
// 绕地球平移（12 天回到原点）。速度场：
//     u = u0 (cosφ cosα + sinφ cosλ sinα)
//     v = -u0 sinλ sinα
// 其中 u0 = 2πa / 12天，α 为平流方向角。
// 用于验证平流方案的误差与收敛阶（L1/L2/Linf 范数）。
// ---------------------------------------------------------------------------
class WilliamsonTC2 : public TestCase {
public:
    // 参数：平流方向角 α（默认 π/4），钟型中心 (λc, φc)，钟型半径 R
    StateReal alpha = 0.7853981633974483;   // π/4
    StateReal lambda_c = -0.25 * 3.141592653589793;  // -3π/4 经度
    StateReal phi_c = 0.0;                   // 纬度
    StateReal bell_radius = 3.0 * 3.141592653589793 / 9.0;  // R = π/3

    void initialize(dynamics::State& state,
                    const grid::VerticalCoordinate& vert) override;
    const char* name() const override { return "WilliamsonTC2"; }

    // 解析解：给定时刻 t（秒）后的标量场（位温/密度示踪剂）
    // 返回钟型场在经纬度 (lon, lat) 处的值（0–1 归一化）。
    StateReal analytic(StateReal lon, StateReal lat, StateReal t) const;

    // 计算相对初始时刻的误差范数（L2），需要网格遍历。
    StateReal error_L2(const dynamics::State& state,
                       const grid::CubedSphereGrid& grid,
                       StateReal t) const;
};

// ---------------------------------------------------------------------------
// Williamson TC5：地形上的罗斯贝波
// ---------------------------------------------------------------------------
// 恒定纬向流跨越孤立高斯山地形，产生定常罗斯贝波解。用于检验地形处理与
// 平衡性。地形高度：
//     h_s(λ, φ) = h0 (1 - r/R)
// 其中 r 为到山中心的大圆距离，R = π/9，h0 = 2000 m。
// ---------------------------------------------------------------------------
class WilliamsonTC5 : public TestCase {
public:
    StateReal h0 = 2000.0;                  // 山高 [m]
    StateReal mountain_lambda = 0.5 * 3.141592653589793;  // 山中心经度 3π/2
    StateReal mountain_phi = 1.0 / 6.0 * 3.141592653589793; // 山中心纬度 π/6
    StateReal mountain_R = 3.141592653589793 / 9.0;  // R = π/9

    void initialize(dynamics::State& state,
                    const grid::VerticalCoordinate& vert) override;
    const char* name() const override { return "WilliamsonTC5"; }

    // 高斯山地形高度
    StateReal surface_height(StateReal lon, StateReal lat) const;
};

// ---------------------------------------------------------------------------
// Jablonowski–Williamson 斜压波（三维）
// ---------------------------------------------------------------------------
// 解析的、静力平衡的斜压波初始条件（Jablonowski & Williamson 2006），
// 含纬向喷流与扰动，用于检验三维动力核心的非静力/斜压响应。
// ---------------------------------------------------------------------------
class JablonowskiWilliamson : public TestCase {
public:
    void initialize(dynamics::State& state,
                    const grid::VerticalCoordinate& vert) override;
    const char* name() const override { return "JablonowskiWilliamson"; }
};

// ---------------------------------------------------------------------------
// Held–Suarez 气候态强迫（三维长期积分）
// ---------------------------------------------------------------------------
// 理想化气候态：温度弛豫到解析平衡态 + 边界层瑞利摩擦。
// 提供温度平衡剖面与摩擦系数，供物理过程插件调用。
// ---------------------------------------------------------------------------
class HeldSuarez : public TestCase {
public:
    void initialize(dynamics::State& state,
                    const grid::VerticalCoordinate& vert) override;
    const char* name() const override { return "HeldSuarez"; }

    // 平衡温度剖面（K），z 为高度 [m]，lat 为纬度 [rad]
    static StateReal equilibrium_temperature(StateReal z, StateReal lat);

    // 瑞利摩擦系数（s⁻¹），z 为高度 [m]
    static StateReal rayleigh_friction(StateReal z);

    // 温度弛豫时间尺度（s）
    static constexpr StateReal relaxation_time = 40.0 * 86400.0;
};

}  // namespace cubed_sph::testcases
