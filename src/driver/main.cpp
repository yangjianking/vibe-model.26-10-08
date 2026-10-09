// =============================================================================
//  主程序入口
// =============================================================================
//  用法：cubed_sphere_nwp [--config <file.yaml>] [--set key=value ...]
//
//  职责（轻量，仅装配与驱动，不含计算内核）：
//   1. 初始化 MPI 并行环境
//   2. 解析配置（YAML + 命令行覆盖）
//   3. 装配 Driver 并运行
//   4. 写 provenance（配置哈希、git hash、编译器版本）
//   5. 终结 MPI 环境
// =============================================================================

#include "cubed_sph/driver/driver.hpp"
#include "cubed_sph/parallel/parallel.hpp"
#include "cubed_sph/common/config.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace {

// 简单的命令行参数解析（--config / --set）
struct Cmdline {
    std::vector<std::string> configs;
    std::vector<std::string> overrides;
};

Cmdline parse_args(int argc, char** argv) {
    Cmdline cmd;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--config" && i + 1 < argc) {
            cmd.configs.push_back(argv[++i]);
        } else if (a == "--set" && i + 1 < argc) {
            cmd.overrides.push_back(argv[++i]);
        } else if (a == "--help" || a == "-h") {
            std::cout << "用法: cubed_sphere_nwp [--config <file.yaml>] "
                         "[--set key=value ...]\n";
        }
    }
    return cmd;
}

}  // namespace

int main(int argc, char** argv) {
    using namespace cubed_sph;

    // 1. 并行环境初始化
    parallel::ParallelContext::init(&argc, &argv);

    // 2. 参数解析与配置加载
    const Cmdline cmd = parse_args(argc, argv);
    std::vector<std::string> paths = cmd.configs;
    if (paths.empty()) {
        paths.push_back("configs/default.yaml");
    }

    try {
        Config cfg = Config::load(paths);
        cfg.apply_cmdline(cmd.overrides);

        // 3. 装配并运行
        Driver driver(cfg);
        driver.initialize();
        driver.run();
        driver.write_output();

        // 4. 输出配置摘要（供 provenance 归档）
        std::cout << "[main] 配置哈希: " << cfg.content_hash() << "\n";
    } catch (const ConfigError& e) {
        std::cerr << "[main] 配置错误: " << e.what() << "\n";
        parallel::ParallelContext::finalize();
        return 1;
    }

    // 5. 并行环境终结
    parallel::ParallelContext::finalize();
    return 0;
}
