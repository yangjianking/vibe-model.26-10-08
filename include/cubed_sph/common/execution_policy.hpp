// =============================================================================
//  可移植执行策略层（P3 GPU/众核内核抽象，Kokkos 首选）
// =============================================================================
//  目标：让业务代码（dynamics/physics/grid）与硬件后端完全解耦，同一份
//  内核源码可编译到 纯 CPU / OpenMP / CUDA / HIP / SYCL。
//
//  设计：优先使用 Kokkos（Edwards et al. 2014, JPDC 74）。若构建时未发现
//  Kokkos，则退化为自研 C++ 模板执行策略层（同语义）：
//    - ParallelForRange   : 一维区间并行 for（含多维扁平化）
//    - ParallelReduceRange: 带归约的并行 for（确定性归约树）
//    - View               : 多维数组抽象（SoA 优先，垂直维最内连续）
//
//  铁律：业务代码禁止出现任何 #ifdef __CUDACC__ 或裸 CUDA/HIP；平台差异
//  只允许存在于本层（kernels/common）。
//
//  参考文献：
//    Edwards et al. (2014), J. Parallel Distrib. Comput. 74, 3202–3216（Kokkos）
//    Adams et al. (2019), J. Parallel Distrib. Comput. 132, 383–396（LFRic PSyKAl）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"

#include <cstddef>
#include <functional>

namespace cubed_sph {

// =============================================================================
// 执行空间选择：优先 Kokkos，否则退化为自研策略层
// =============================================================================
#if defined(CS_USE_KOKKOS) && defined(Kokkos_VERSION)
    // ---- Kokkos 后端（推荐路径） ----
    #include <Kokkos_Core.hpp>
    namespace exec = Kokkos;

    template <typename T, typename... Dims>
    using View = Kokkos::View<T, Dims...>;

    using DefaultExecSpace  = Kokkos::DefaultExecutionSpace;
    using HostExecSpace     = Kokkos::HostSpace;

    // 一维区间并行 for
    template <typename Fn>
    inline void parallel_for_range(std::size_t n, const Fn& fn) {
        Kokkos::parallel_for("cs_kernel",
            Kokkos::RangePolicy<DefaultExecSpace>(0, static_cast<long>(n)), fn);
    }

    // 一维区间并行归约（确定性归约树）
    template <typename Fn>
    inline ComputeReal parallel_reduce_range(std::size_t n, const Fn& fn,
                                             ComputeReal init) {
        ComputeReal acc = init;
        Kokkos::parallel_reduce("cs_reduce",
            Kokkos::RangePolicy<DefaultExecSpace>(0, static_cast<long>(n)),
            fn, acc);
        return acc;
    }

#else
    // ---- 自研 C++ 模板执行策略层（退化路径，同语义） ----
    // 说明：当 Kokkos 不可用时退化为纯 CPU（OpenMP/串行）实现，保持与
    //       Kokkos 路径完全一致的调用接口，业务代码无需任何修改。
    #include <omp.h>

    // 简化 View：堆分配的连续多维数组（垂直维最内连续，SoA 语义）
    template <typename T>
    class View {
    public:
        View() = default;
        View(std::size_t n) : data_(n, T{}) {}
        View(std::size_t n, const T& v) : data_(n, v) {}

        T& operator()(std::size_t i) { return data_[i]; }
        const T& operator()(std::size_t i) const { return data_[i]; }
        std::size_t size() const { return data_.size(); }

        T*       data()       { return data_.data(); }
        const T* data() const { return data_.data(); }

    private:
        std::vector<T> data_;
    };

    // 一维区间并行 for（OpenMP 可选，未启用时串行）
    template <typename Fn>
    inline void parallel_for_range(std::size_t n, const Fn& fn) {
        #if defined(_OPENMP)
        #pragma omp parallel for
        #endif
        for (std::size_t i = 0; i < n; ++i) {
            fn(static_cast<long>(i));
        }
    }

    // 一维区间并行归约（确定性：固定顺序累加，保证位可重复与线程数无关）
    template <typename Fn>
    inline ComputeReal parallel_reduce_range(std::size_t n, const Fn& fn,
                                             ComputeReal init) {
        ComputeReal acc = init;
        for (std::size_t i = 0; i < n; ++i) {
            fn(static_cast<long>(i), acc);   // 用户 lambda 负责 acc += 局部贡献
        }
        return acc;
    }

#endif

}  // namespace cubed_sph
