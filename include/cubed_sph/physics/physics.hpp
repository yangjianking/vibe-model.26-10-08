// =============================================================================
//  物理参数化接口（Physics Parameterization Interface）
// =============================================================================
//  物理过程（辐射、边界层、微物理、对流、重力波拖曳、地表过程）通过统一
//  插件接口接入动力核心。每个参数化方案实现 PhysicalProcess 接口，向模式
//  提供：
//    - tendency(state, t)    ：单位时间倾向（加到位温、动量、水汽等）
//    - name()                ：方案名（用于日志与 provenance）
//
//  设计原则（可扩展、可测试）：
//    - 插件式注册：新方案只需实现接口并在工厂注册，无需改动动力核心；
//    - 倾向累加：多个方案按配置顺序依次调用，倾向线性叠加；
//    - 时间分离：物理倾向通常用较长时间步（与动力快波解耦，operator splitting）。
//
//  首批方案（占位实现，接口完整）：
//    - DummyRadiation      ：辐射参数化占位
//    - DummyBoundaryLayer  ：边界层参数化占位
//    - DummyMicrophysics   ：微物理参数化占位
//    - DummyGravityWaveDrag：重力波拖曳占位
//
//  参考文献：
//    Held & Suarez (1994), BAMS 75, 1825–1830（简化物理，气候态验证）
//    Emanuel (1994), Atmospheric Convection（对流参数化理论）
// =============================================================================

#pragma once

#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/common/config.hpp"

#include <string>
#include <vector>
#include <memory>

namespace cubed_sph::physics {

// ---------------------------------------------------------------------------
// 物理过程抽象接口
// ---------------------------------------------------------------------------
class PhysicalProcess {
public:
    virtual ~PhysicalProcess() = default;

    // 计算物理倾向并累加到 tendency（dt 为物理时间步长）
    virtual void tendency(dynamics::State& state, StateReal dt) = 0;

    // 方案名（用于日志、provenance、配置校验）
    virtual const char* name() const = 0;
};

// ---------------------------------------------------------------------------
// 物理包管理器：按配置注册方案并统一调度
// ---------------------------------------------------------------------------
class PhysicsPackage {
public:
    PhysicsPackage() = default;

    // 注册一个参数化方案（工厂式，供 config 驱动）
    void add(std::unique_ptr<PhysicalProcess> process);

    // 依次调用所有已注册方案的倾向
    void apply(dynamics::State& state, StateReal dt);

    // 已注册方案名列表（用于日志）
    std::vector<std::string> active_processes() const;

private:
    std::vector<std::unique_ptr<PhysicalProcess>> processes_;
};

// ---------------------------------------------------------------------------
// 占位方案（接口完整、内部无物理，供框架联调；后续替换为真实参数化）
// ---------------------------------------------------------------------------
class DummyRadiation : public PhysicalProcess {
public:
    void tendency(dynamics::State& state, StateReal dt) override;
    const char* name() const override { return "DummyRadiation"; }
};

class DummyBoundaryLayer : public PhysicalProcess {
public:
    void tendency(dynamics::State& state, StateReal dt) override;
    const char* name() const override { return "DummyBoundaryLayer"; }
};

class DummyMicrophysics : public PhysicalProcess {
public:
    void tendency(dynamics::State& state, StateReal dt) override;
    const char* name() const override { return "DummyMicrophysics"; }
};

}  // namespace cubed_sph::physics
