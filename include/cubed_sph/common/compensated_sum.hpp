// =============================================================================
//  补偿求和（P3 混合精度归约，禁朴素求和）
// =============================================================================
//  全局守恒诊断、能量/质量积分等归约必须使用补偿求和，避免朴素累加的
//  灾难性抵消误差。提供 Kahan 与 Neumaier 两种算法。
//
//  参考文献：
//    Kahan, W. (1965), Commun. ACM 8, 40（Kahan 补偿求和）
//    Higham, N. J. (2002), Accuracy and Stability of Numerical Algorithms,
//      2nd ed., SIAM（Neumaier 变体，§4.3）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"

namespace cubed_sph {

// ---------------------------------------------------------------------------
// Kahan 补偿求和累加器
//  - 误差项 c 补偿低位舍入，使累计误差接近 O(ε) 而非 O(n·ε)。
// ---------------------------------------------------------------------------
class KahanSum {
public:
    void add(ComputeReal x) {
        const ComputeReal y = x - c_;
        const ComputeReal t = sum_ + y;
        c_ = (t - sum_) - y;   // 恢复丢失的低位
        sum_ = t;
    }
    ComputeReal value() const { return sum_; }
    void reset() { sum_ = 0.0; c_ = 0.0; }

private:
    ComputeReal sum_ = 0.0;
    ComputeReal c_   = 0.0;   // 补偿项
};

// ---------------------------------------------------------------------------
// Neumaier 补偿求和累加器
//  - 相比 Kahan，当 |sum| ≥ |x| 时能额外补偿，对单精度累加更稳健。
// ---------------------------------------------------------------------------
class NeumaierSum {
public:
    void add(ComputeReal x) {
        const ComputeReal t = sum_ + x;
        if (std::abs(sum_) >= std::abs(x)) {
            c_ += (sum_ - t) + x;   // 低阶项进入补偿
        } else {
            c_ += (x - t) + sum_;
        }
        sum_ = t;
    }
    ComputeReal value() const { return sum_ + c_; }
    void reset() { sum_ = 0.0; c_ = 0.0; }

private:
    ComputeReal sum_ = 0.0;
    ComputeReal c_   = 0.0;
};

}  // namespace cubed_sph
