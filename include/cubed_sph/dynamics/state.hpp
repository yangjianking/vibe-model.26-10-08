// =============================================================================
//  模式状态（State）—— 预后变量与诊断变量的容器
// =============================================================================
//  非静力可压缩 Euler 方程的预后变量（提示词 P2 第 1 条）：
//    ρ   ：干空气质量密度（或总密度，视水汽处理而定）
//    ρu  ：ξ 方向动量密度（协变分量存储，见 grid/covariant.hpp）
//    ρv  ：η 方向动量密度
//    ρw  ：垂直动量密度
//    ρθ  ：位温密度（θ 为位温；或用 ρπ 形式的 Exner 压力，见 thermo.hpp）
//
//  采用"密度×速度"形式的动机：连续方程与动量方程天然写成通量守恒形式，
//  便于散度/压力梯度严格配对（Arakawa–Lamb 1981），保证离散总质量/动量
//  守恒。
//
//  存储布局：与网格一致（panel × n × n），垂直维 k 最内连续（SoA 优先，
//  利于 cache 与 GPU coalescing，见提示词 P3 第 3 条）。
//
//  参考文献：
//    Harris et al. (2021), GFDL TM GFDL2021001（FV³ 变量组织）
//    Wood et al. (2014), QJRMS 140（ENDGame 变量与守恒）
// =============================================================================

#pragma once

#include "cubed_sph/grid/cubed_sphere_grid.hpp"
#include "cubed_sph/common/types.hpp"

#include <vector>

namespace cubed_sph::dynamics {

// ---------------------------------------------------------------------------
// 状态字段：以扁平 vector 存储，索引 = panel*n*n*nlev + k*n*n + j*n + i
// ---------------------------------------------------------------------------
class State {
public:
    // 构造函数：根据网格与垂直层数分配
    State(const grid::CubedSphereGrid& grid, int nlev);

    // 拷贝赋值：仅深拷贝字段数据（grid_ 为引用，保持指向同一网格）。
    // 注意：State 含引用成员，编译器生成的 operator= 会被删除，
    // 故显式实现。要求两侧绑定同一网格（由调用方保证）。
    State& operator=(const State& other);

    // 维度访问
    int nlev() const { return nlev_; }
    const grid::CubedSphereGrid& grid() const { return grid_; }
    IIndex ncell() const { return static_cast<IIndex>(nlev_) * ncell2d_; }

    // 单个字段的扁平索引（panel, i, j, k）
    IIndex index(int panel, int i, int j, int k) const;

    // 预后变量字段（协变动量分量）
    std::vector<StateReal> rho;    // 密度
    std::vector<StateReal> rho_u;  // ρu（ξ 协变）
    std::vector<StateReal> rho_v;  // ρv（η 协变）
    std::vector<StateReal> rho_w;  // ρw（垂直）
    std::vector<StateReal> rho_theta;  // ρθ（位温密度）

    // 诊断/辅助字段
    std::vector<StateReal> exner;   // Exner 压力 π = (p/p0)^κ（可选，见 thermo）
    std::vector<StateReal> theta;   // 位温 θ
    std::vector<StateReal> geopotential;  // 重力位 Φ = gz

    // 派生速度（由动量/密度还原，供诊断输出）
    void diagnose_velocities();

private:
    const grid::CubedSphereGrid& grid_;
    int nlev_ = 1;
    IIndex ncell2d_ = 0;  // 单层水平单元数（含 halo，6*n*n）
};

}  // namespace cubed_sph::dynamics
