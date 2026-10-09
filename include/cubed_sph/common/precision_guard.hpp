// =============================================================================
//  混合精度抽象：精度卫士与随机舍入（P3 第 4 条）
// =============================================================================
//  三级类型别名 Real / ComputeReal / StateReal 已在 types.hpp 定义。本模块
//  提供混合精度运行时的两个关键辅助：
//
//   1. 精度卫士（precision guard）：除零/上下溢保护，避免 fp32 降精度导致
//      的 NaN/Inf 污染。参考 Váňa et al. (2017) 的经验。
//   2. 随机舍入（stochastic rounding）：可选开关，用于精度敏感性诊断
//      （Paxton et al. 2022）。
//
//  误差预算策略（默认，见 docs/numerics/mixed_precision.md）：
//    - 网格几何量、椭圆求解器、全局守恒诊断、状态累积量 → fp64（StateReal）
//    - 平流通量、物理倾向、扩散等主体计算 → fp32（Real）
//    - 归约用 Kahan/Neumaier 补偿求和（compensated_sum.hpp）
//
//  参考文献：
//    Váňa et al. (2017), Mon. Wea. Rev. 145, 4957–4970（单精度评估）
//    Düben & Palmer (2014), Mon. Wea. Rev. 142, 3809–3826（不精确硬件）
//    Paxton et al. (2022) [TO-VERIFY]（随机舍入）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"

#include <cmath>
#include <limits>
#include <cstdint>

namespace cubed_sph::common {

// ---------------------------------------------------------------------------
// 精度卫士：除零保护
// ---------------------------------------------------------------------------
// 分母加一个极小正数，避免 fp32 下除以 0 或极小值产生 Inf。
// 阈值取 fp32 最小正规数的量级（~1e-38）的若干倍，兼顾精度与稳健。
inline StateReal guard_divisor(StateReal denom,
                               StateReal floor = 1.0e-20) {
    return (std::abs(denom) < floor)
               ? (denom >= 0.0 ? floor : -floor)
               : denom;
}

// ---------------------------------------------------------------------------
// 精度卫士：有限性检查（NaN/Inf 检测，用于降精度点诊断）
// ---------------------------------------------------------------------------
inline bool is_finite_value(StateReal x) {
    return std::isfinite(static_cast<double>(x));
}

// ---------------------------------------------------------------------------
// 随机舍入（stochastic rounding，可选开关）
// ---------------------------------------------------------------------------
// 将 fp64 值舍入到 fp32 时，以随机方式决定进位（期望值等于精确值），
// 用于消除系统性舍入偏差，评估精度敏感性（Paxton et al. 2022）。
//
// 注意：随机舍入会破坏位可重复性，仅用于诊断，生产运行默认关闭。
class StochasticRounder {
public:
    // 用给定的整数种子初始化随机数生成器（确定性，便于复现诊断）
    explicit StochasticRounder(std::uint64_t seed = 12345) : state_(seed) {}

    // 将 fp64 舍入到 fp32（随机舍入），期望值 = x
    float round_to_f32(double x) {
        const float lo = static_cast<float>(x);
        const float hi = std::nextafter(lo,
            x >= 0.0 ? std::numeric_limits<float>::infinity()
                     : -std::numeric_limits<float>::infinity());
        // 概率 = (x - lo) / (hi - lo)
        const double frac = (x - static_cast<double>(lo)) /
                            (static_cast<double>(hi) - static_cast<double>(lo) + 1e-300);
        const double u = next_uniform();
        return (u < frac) ? hi : lo;
    }

private:
    std::uint64_t state_;

    // xorshift64 均匀随机数生成器（确定性）
    double next_uniform() {
        state_ ^= state_ << 13;
        state_ ^= state_ >> 7;
        state_ ^= state_ << 17;
        const std::uint64_t r = state_;
        return static_cast<double>(r >> 11) * (1.0 / 9007199254740992.0);
    }
};

}  // namespace cubed_sph::common
