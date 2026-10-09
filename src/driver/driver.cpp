// =============================================================================
//  模式驱动实现
// =============================================================================

#include "cubed_sph/driver/driver.hpp"

#include <iostream>

namespace cubed_sph {

Driver::Driver(const Config& cfg) : cfg_(cfg) {
    parse_grid_params();
    parse_time_params();
}

void Driver::parse_grid_params() {
    ncells_ = cfg_.get_int("grid.ncells", 48);
    nhalo_  = cfg_.get_int("grid.nhalo", 3);
    nlev_   = cfg_.get_int("vertical.nlev", 60);
}

void Driver::parse_time_params() {
    dt_ = cfg_.get_double("timeint.dt", 300.0);
    nsteps_ = cfg_.get_int("run.nsteps", 10);
}

void Driver::initialize() {
    // 1. 构建立方球网格（含度量缓存）
    grid::GridParams gp;
    gp.ncells = ncells_;
    gp.nhalo = nhalo_;
    grid_ = std::make_unique<grid::CubedSphereGrid>(gp);

    // 2. 构建垂直坐标（地形跟随 / 纯高度，由配置选择）
    grid::VerticalParams vp;
    vp.nlev = nlev_;
    vp.z_top = cfg_.get_double("vertical.model_top", 30000.0);
    vp.s = cfg_.get_double("vertical.smoothing_s", 5.0);
    const std::string vtype = cfg_.get_string("vertical.type", "terrain_following");
    vp.type = (vtype == "pure_height") ? grid::VerticalCoordType::PureHeight
                                       : grid::VerticalCoordType::TerrainFollowing;
    vert_ = std::make_unique<grid::VerticalCoordinate>(vp);

    // 3. 初始化状态
    state_ = std::make_unique<dynamics::State>(*grid_, nlev_);

    // 4. 构建动力核心
    core_ = std::make_unique<dynamics::DynamicsCore>(*grid_, *vert_, nlev_, cfg_);

    // 5. 构建时间积分器（半隐式 or 显式，由配置选择）
    timeint::TimeStepperParams tsp;
    tsp.dt = dt_;
    tsp.off_centering = cfg_.get_double("timeint.off_centering", 0.5);
    tsp.use_semi_implicit = cfg_.get_bool("timeint.semi_implicit", true);
    if (tsp.use_semi_implicit) {
        integrator_ = std::make_unique<timeint::SemiImplicitIntegrator>(
            *core_, *vert_, tsp);
    } else {
        integrator_ = std::make_unique<timeint::ExplicitRK3Integrator>(*core_, tsp);
    }

    // 6. 构建物理包（由配置注册方案）
    build_physics();

    std::cout << "[Driver] 初始化完成：ncells=" << ncells_
              << ", nhalo=" << nhalo_ << ", nlev=" << nlev_
              << ", dt=" << dt_ << "s, steps=" << nsteps_ << "\n";
    std::cout << "[Driver] 垂直坐标："
              << (vp.type == grid::VerticalCoordType::TerrainFollowing
                      ? "地形跟随" : "纯高度")
              << ", 模式顶=" << vp.z_top << "m\n";
    std::cout << "[Driver] 时间积分器：" << integrator_->name() << "\n";
}

void Driver::build_physics() {
    physics_ = std::make_unique<physics::PhysicsPackage>();
    // 由配置开关注册占位方案（后续替换为真实参数化）
    if (cfg_.get_bool("physics.radiation", false)) {
        physics_->add(std::make_unique<physics::DummyRadiation>());
    }
    if (cfg_.get_bool("physics.boundary_layer", false)) {
        physics_->add(std::make_unique<physics::DummyBoundaryLayer>());
    }
    if (cfg_.get_bool("physics.microphysics", false)) {
        physics_->add(std::make_unique<physics::DummyMicrophysics>());
    }
}

void Driver::run() {
    std::cout << "[Driver] 开始时间推进（" << nsteps_ << " 步）\n";
    for (int step = 0; step < nsteps_; ++step) {
        integrator_->step(*state_);
        // 物理过程（operator splitting，此处占位调度）
        physics_->apply(*state_, dt_);
    }
    std::cout << "[Driver] 时间推进完成\n";
}

void Driver::write_output() {
    // TODO(driver): 装配 io::NetCDFWriter 输出预报场。
    std::cout << "[Driver] 输出接口已预留（见 src/io/）\n";
}

}  // namespace cubed_sph
