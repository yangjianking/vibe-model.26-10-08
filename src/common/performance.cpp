// =============================================================================
//  性能基准与计时实现
// =============================================================================

#include "cubed_sph/common/performance.hpp"

#include <sstream>
#include <algorithm>

namespace cubed_sph::common {

// ---------------------------------------------------------------------------
// KernelTimer
// ---------------------------------------------------------------------------
void KernelTimer::record(const std::string& kernel_name, double sec) {
    records_.emplace_back(kernel_name, sec);
    total_ += sec;
}

std::string KernelTimer::report() const {
    std::ostringstream oss;
    for (const auto& [name, sec] : records_) {
        oss << name << ": " << sec << " s\n";
    }
    oss << "total: " << total_ << " s\n";
    return oss.str();
}

void KernelTimer::reset() {
    total_ = 0.0;
    records_.clear();
}

// ---------------------------------------------------------------------------
// RooflineBench
// ---------------------------------------------------------------------------
double RooflineBench::stream_bandwidth_gbps(std::size_t n, int n_repeat) {
    // STREAM Triad：a[i] = b[i] + s * c[i]
    // 每个元素：读 2（b,c）+ 写 1（a）= 3 次内存访问 × 8 字节 = 24 字节
    std::vector<StateReal> a(n, 1.0), b(n, 2.0), c(n, 0.5);
    const StateReal s = 3.0;

    // 预热（消除首次访问的冷缓存）
    for (std::size_t i = 0; i < n; ++i) a[i] = b[i] + s * c[i];

    Timer t;
    t.start();
    for (int r = 0; r < n_repeat; ++r) {
        for (std::size_t i = 0; i < n; ++i) {
            a[i] = b[i] + s * c[i];
        }
    }
    const double sec = t.elapsed_sec();

    const double bytes = 3.0 * static_cast<double>(n) * sizeof(StateReal) *
                         static_cast<double>(n_repeat);
    const double gb = bytes / 1.0e9;
    return gb / sec;
}

double RooflineBench::peak_flops(std::size_t n, int n_repeat) {
    // 简化的 FMA 吞吐测试：acc += a[i] * b[i]（1 FMA = 2 FLOP）
    std::vector<StateReal> a(n, 1.5), b(n, 2.5);
    StateReal acc = 0.0;

    // 预热
    for (std::size_t i = 0; i < n; ++i) acc += a[i] * b[i];

    Timer t;
    t.start();
    for (int r = 0; r < n_repeat; ++r) {
        StateReal local = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            local += a[i] * b[i];  // 1 FMA = 2 FLOP
        }
        acc += local;
    }
    const double sec = t.elapsed_sec();

    const double flops = 2.0 * static_cast<double>(n) *
                         static_cast<double>(n_repeat);
    (void)acc;  // 防止优化掉
    return flops / sec / 1.0e9;  // GFLOPS
}

double RooflineBench::arithmetic_intensity(double flops, double bytes) {
    return flops / bytes;  // FLOP/Byte
}

}  // namespace cubed_sph::common
