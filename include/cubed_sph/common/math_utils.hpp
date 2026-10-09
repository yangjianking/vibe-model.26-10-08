// =============================================================================
//  数学工具（球面几何、向量旋转、平滑函数）
// =============================================================================
//  集中放置跨模块共享的数学原语，避免各处重复实现导致的不一致。
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"

#include <array>
#include <cmath>

namespace cubed_sph::math {

// ---------------------------------------------------------------------------
// 三维单位向量（笛卡尔分量，用于面板间向量旋转）
// ---------------------------------------------------------------------------
struct Vec3 {
    ComputeReal x = 0.0;
    ComputeReal y = 0.0;
    ComputeReal z = 0.0;
};

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x};
}

inline ComputeReal dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 normalize(const Vec3& v) {
    const ComputeReal n = std::sqrt(dot(v, v));
    return Vec3{v.x / n, v.y / n, v.z / n};
}

// ---------------------------------------------------------------------------
// 平滑阶跃 / 权重函数（嵌套边界松弛、海绵层、扩散系数标定）
// ---------------------------------------------------------------------------
// 参考：Davies (1976, QJRMS 102) 侧边界松弛；Klemp et al. (2008, MWR 136)
// 模顶 Rayleigh 海绵层。平滑权重的通用形式为 [0,1] 上的单调 C¹ 函数。

// 三次多项式平滑（端点为 0/1 且一阶导为 0）
inline ComputeReal smooth_cubic(ComputeReal t) {
    // t ∈ [0,1]，返回 3t² - 2t³，C¹ 连续
    const ComputeReal tt = t * t;
    return tt * (3.0 - 2.0 * t);
}

// 双曲正切平滑（可配过渡宽度）
inline ComputeReal smooth_tanh(ComputeReal t, ComputeReal width) {
    // width 控制过渡区宽度，越小越陡
    return 0.5 * (1.0 + std::tanh((t - 0.5) / (width + 1e-12)));
}

// ---------------------------------------------------------------------------
// 球面距离（两单位向量夹角，用于嵌套松弛系数随距离衰减）
// ---------------------------------------------------------------------------
inline ComputeReal spherical_distance(const Vec3& a, const Vec3& b) {
    const ComputeReal c = dot(a, b);
    const ComputeReal clamped = std::max(-1.0, std::min(1.0, c));
    return std::acos(clamped);
}

}  // namespace cubed_sph::math
