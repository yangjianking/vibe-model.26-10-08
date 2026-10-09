// =============================================================================
//  SoA 内存布局（Structure-of-Arrays，垂直维最内连续）
// =============================================================================
//  提示词 P3 第 3 条：SoA 优先；垂直维为最内连续维以利 cache 与 GPU coalescing；
//  提供 AoS/SoA/SoA-tile 编译期切换便于基准对比。
//
//  布局约定：
//    - 水平索引 h ∈ [0, N_h)（N_h = 6 × n × n，含 halo）
//    - 垂直索引 k ∈ [0, nlev)
//    - SoA 布局：field[h * nlev + k]，垂直维连续（k 相邻 → 内存相邻）
//    这保证同一水平格点上的垂直列在内存中连续，利于：
//      * 垂直隐式求解（三对角/线松弛）的 cache 局部性；
//      * GPU 线程按垂直维 coalesced 访问。
//
//  与现有 State（AoS：field[panel][k][j][i]）的对比：
//    - State 当前用 AoS（k 为中间维）；本模块提供 SoA 容器与互转，
//      供性能基准对比（提示词 P3 验收的 thread-count invariance 与带宽测试）。
//
//  参考文献：
//    Edwards et al. (2014), JPDC 74（Kokkos View 布局）
//    Harris et al. (2021), GFDL TM GFDL2021001 §2（FV³ 内存布局）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"
#include "cubed_sph/common/execution_policy.hpp"

#include <vector>

namespace cubed_sph::common {

// ---------------------------------------------------------------------------
// SoA 字段容器：一个变量（如 ρ、ρu）的全部格点，垂直维最内连续。
// ---------------------------------------------------------------------------
// 布局：field[h * nlev + k]，其中 h 为水平扁平索引，k 为垂直索引。
// 相比 AoS（field[h + k * Nh]），SoA 让垂直列连续。
class SoAField {
public:
    SoAField() = default;
    SoAField(IIndex n_horizontal, int nlev)
        : nh_(n_horizontal), nlev_(nlev),
          data_(static_cast<size_t>(n_horizontal * nlev), 0.0) {}

    // 访问 field(h, k)
    StateReal& operator()(IIndex h, int k) {
        return data_[static_cast<size_t>(h * nlev_ + k)];
    }
    const StateReal& operator()(IIndex h, int k) const {
        return data_[static_cast<size_t>(h * nlev_ + k)];
    }

    IIndex n_horizontal() const { return nh_; }
    int nlev() const { return nlev_; }
    StateReal* data() { return data_.data(); }
    const StateReal* data() const { return data_.data(); }

    void fill(StateReal v) { std::fill(data_.begin(), data_.end(), v); }

private:
    IIndex nh_ = 0;
    int nlev_ = 0;
    std::vector<StateReal> data_;
};

// ---------------------------------------------------------------------------
// AoS ↔ SoA 互转（供性能基准与诊断）
// ---------------------------------------------------------------------------
// 从 AoS 布局（State 当前布局，索引 = panel*n*n*nlev + k*n*n + j*n + i）
// 转为 SoA（索引 = h*nlev + k，h = panel*n*n + j*n + i）。
// ---------------------------------------------------------------------------
inline void aos_to_soa(const std::vector<StateReal>& aos,
                       int panel_n, int nlev,
                       SoAField& soa) {
    const int n = panel_n;
    for (int p = 0; p < 6; ++p) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                const IIndex h = static_cast<IIndex>(p) * n * n +
                                 static_cast<IIndex>(j) * n + i;
                for (int k = 0; k < nlev; ++k) {
                    const IIndex a = static_cast<IIndex>(p) * n * n * nlev +
                                     static_cast<IIndex>(k) * n * n +
                                     static_cast<IIndex>(j) * n + i;
                    soa(h, k) = aos[a];
                }
            }
        }
    }
}

inline void soa_to_aos(const SoAField& soa,
                       int panel_n, int nlev,
                       std::vector<StateReal>& aos) {
    const int n = panel_n;
    for (int p = 0; p < 6; ++p) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                const IIndex h = static_cast<IIndex>(p) * n * n +
                                 static_cast<IIndex>(j) * n + i;
                for (int k = 0; k < nlev; ++k) {
                    const IIndex a = static_cast<IIndex>(p) * n * n * nlev +
                                     static_cast<IIndex>(k) * n * n +
                                     static_cast<IIndex>(j) * n + i;
                    aos[a] = soa(h, k);
                }
            }
        }
    }
}

}  // namespace cubed_sph::common
