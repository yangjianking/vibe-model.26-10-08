// =============================================================================
//  性能基准与计时（P3 第 5 条：roofline 微基准 + 内核计时上报）
// =============================================================================
//  提供：
//    - Timer：高精度计时器（内核级耗时上报，每步输出各内核耗时到日志）
//    - RooflineBench：STREAM 型带宽 / DGEMM 型算力微基准，评估硬件 roofline
//    - 确定性归约校验接口（thread-count invariance）
//
//  参考文献：
//    McCalpin, J. (1995). STREAM: Sustainable Memory Bandwidth（带宽基准）
//    Williams, S., Waterman, A., Patterson, D. (2009), Commun. ACM 52, 65–76
//      （roofline 模型）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"

#include <chrono>
#include <string>
#include <vector>

namespace cubed_sph::common {

// ---------------------------------------------------------------------------
// 高精度计时器（单调时钟）
// ---------------------------------------------------------------------------
class Timer {
public:
    void start() { t0_ = std::chrono::steady_clock::now(); }
    // 返回自 start 起的耗时（秒）
    double elapsed_sec() const {
        const auto t1 = std::chrono::steady_clock::now();
        return std::chrono::duration<double>(t1 - t0_).count();
    }

private:
    std::chrono::steady_clock::time_point t0_;
};

// ---------------------------------------------------------------------------
// 内核计时上报（每步累计各内核耗时，输出 YAML 友好日志）
// ---------------------------------------------------------------------------
class KernelTimer {
public:
    // 记录一次内核耗时
    void record(const std::string& kernel_name, double sec);

    // 累计总耗时
    double total_sec() const { return total_; }

    // 输出报告（"kernel: sec" 每行，供 YAML 日志）
    std::string report() const;

    void reset();

private:
    double total_ = 0.0;
    std::vector<std::pair<std::string, double>> records_;
};

// ---------------------------------------------------------------------------
// Roofline 微基准：带宽（STREAM 型）与算力（DGEMM 型）
// ---------------------------------------------------------------------------
class RooflineBench {
public:
    // 内存带宽测试：对长度为 n 的数组做 STREAM Triad 型操作
    //   a[i] = b[i] + s * c[i]，返回实测带宽 [GB/s]
    static double stream_bandwidth_gbps(std::size_t n, int n_repeat = 20);

    // 计算峰值（简化的 FMA 吞吐测试），返回 [GFLOPS]
    static double peak_flops(std::size_t n, int n_repeat = 20);

    // 计算运算强度（arithmetic intensity）与 roofline 拐点
    static double arithmetic_intensity(double flops, double bytes);
};

}  // namespace cubed_sph::common
