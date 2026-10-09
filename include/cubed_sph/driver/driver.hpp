// =============================================================================
//  模式驱动（Driver）—— 装配各模块并运行主循环
// =============================================================================
//  职责：解析配置 → 构建网格 → 初始化状态 → 装配动力/物理/时间积分/I/O →
//  运行时间推进主循环 → 定期输出与重启。
//
//  这是把 grid/dynamics/timeint/physics/io/parallel 各模块装配成完整可执行
//  程序的编排层（不包含任何计算内核，符合 P0 职责分离）。
// =============================================================================

#pragma once

#include "cubed_sph/common/config.hpp"
#include "cubed_sph/grid/cubed_sphere_grid.hpp"
#include "cubed_sph/grid/vertical_coordinate.hpp"
#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/dynamics/equations.hpp"
#include "cubed_sph/timeint/time_integrator.hpp"
#include "cubed_sph/physics/physics.hpp"

#include <memory>

namespace cubed_sph {

// ---------------------------------------------------------------------------
// 模式驱动
// ---------------------------------------------------------------------------
class Driver {
public:
    // 从配置构建驱动（不启动，仅装配）
    explicit Driver(const Config& cfg);

    // 初始化各模块（分配内存、构建网格度量、装配时间积分与物理包）
    void initialize();

    // 运行主循环（nsteps 步，或由配置指定）
    void run();

    // 输出当前状态到文件
    void write_output();

private:
    Config cfg_;

    // 网格参数（由配置解析）
    int ncells_ = 48;
    int nhalo_ = 3;
    int nlev_ = 60;
    int nsteps_ = 0;
    StateReal dt_ = 300.0;

    // 模块实例
    std::unique_ptr<grid::CubedSphereGrid> grid_;
    std::unique_ptr<grid::VerticalCoordinate> vert_;
    std::unique_ptr<dynamics::State> state_;
    std::unique_ptr<dynamics::DynamicsCore> core_;
    std::unique_ptr<timeint::TimeIntegrator> integrator_;
    std::unique_ptr<physics::PhysicsPackage> physics_;

    void parse_grid_params();
    void parse_time_params();
    void build_physics();
};

}  // namespace cubed_sph
