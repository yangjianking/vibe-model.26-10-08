// =============================================================================
//  物理参数化实现（占位方案 + 管理器）
// =============================================================================

#include "cubed_sph/physics/physics.hpp"

namespace cubed_sph::physics {

void PhysicsPackage::add(std::unique_ptr<PhysicalProcess> process) {
    processes_.push_back(std::move(process));
}

void PhysicsPackage::apply(dynamics::State& state, StateReal dt) {
    for (auto& p : processes_) {
        p->tendency(state, dt);
    }
}

std::vector<std::string> PhysicsPackage::active_processes() const {
    std::vector<std::string> names;
    names.reserve(processes_.size());
    for (const auto& p : processes_) names.push_back(p->name());
    return names;
}

// ---------------------------------------------------------------------------
// 占位方案：不做任何物理，仅保持接口连通。
// 后续替换为真实参数化时，只需实现 tendency() 内部逻辑。
// ---------------------------------------------------------------------------
void DummyRadiation::tendency(dynamics::State& state, StateReal dt) {
    // TODO(physics): 实现长短波辐射参数化（如 RRTMG 接口）。
    // 本阶段占位：辐射加热率为 0。
    (void)state; (void)dt;
}

void DummyBoundaryLayer::tendency(dynamics::State& state, StateReal dt) {
    // TODO(physics): 实现边界层参数化（如 MYNN/YSU）。
    (void)state; (void)dt;
}

void DummyMicrophysics::tendency(dynamics::State& state, StateReal dt) {
    // TODO(physics): 实现微物理参数化（如 Thompson/WSM6）。
    (void)state; (void)dt;
}

}  // namespace cubed_sph::physics
