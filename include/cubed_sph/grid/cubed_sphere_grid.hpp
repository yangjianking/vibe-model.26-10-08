// =============================================================================
//  立方球网格（CubedSphereGrid）—— 水平网格主类
// =============================================================================
//  封装 6 个面板的网格拓扑、度量缓存、Arakawa C/D 交错索引，以及跨面板
//  halo 交换与向量旋转。这是动力学核心与 I/O 子系统共享的网格数据源。
//
//  网格规模参数（分辨率约定）：
//    - ncells：每个面板一个方向上的单元数（不含 halo）。全球水平单元总数
//      ≈ 6 × ncells²。常见分辨率对应：Cs48 → ncells=48；Cs96 → ncells=96。
//    - nhalo ：halo 层数（默认 3，覆盖平流 + 隐式扩散 + 嵌套边界缓冲）。
//
//  交错约定（Arakawa C-grid，水平）：
//    - 标量（ρ, θ, π）    ：单元中心 (i, j)
//    - u（ξ 方向风速）    ：单元东面 (i+1/2, j)
//    - v（η 方向风速）    ：单元北面 (i, j+1/2)
//    - 涡度/散度          ：单元顶点 (i+1/2, j+1/2)（D-grid 可选）
//
//  参考文献：
//    Ronchi et al. (1996), J. Comput. Phys. 124
//    Putman & Lin (2007), J. Comput. Phys. 227, 55–78
//    Harris et al. (2021), GFDL TM GFDL2021001
// =============================================================================

#pragma once

#include "cubed_sph/grid/panel.hpp"
#include "cubed_sph/grid/metric.hpp"
#include "cubed_sph/grid/panel_rotation.hpp"
#include "cubed_sph/common/types.hpp"
#include "cubed_sph/common/execution_policy.hpp"

#include <vector>

namespace cubed_sph::grid {

// ---------------------------------------------------------------------------
// 网格构建参数
// ---------------------------------------------------------------------------
struct GridParams {
    int ncells = 48;     // 每面板单方向单元数（不含 halo）
    int nhalo  = 3;      // halo 层数（默认 3）
    StateReal radius = 6371229.0;  // 地球半径（米），默认与 GFDL FV³ 一致
};

// ---------------------------------------------------------------------------
// 立方球网格
// ---------------------------------------------------------------------------
class CubedSphereGrid {
public:
    explicit CubedSphereGrid(const GridParams& params);

    // ---- 基本维度 ----
    int ncells() const { return ncells_; }
    int nhalo()  const { return nhalo_; }
    int num_panels() const { return kNumPanels; }
    StateReal radius() const { return radius_; }

    // 单面板逻辑尺寸（含 halo）：n = ncells + 2*nhalo
    int panel_n() const { return ncells_ + 2 * nhalo_; }

    // ---- 度量访问 ----
    // 返回面板 p、逻辑索引 (i,j) 处的度量（i,j 为含 halo 的逻辑索引）
    const MetricPoint& metric(int panel, int i, int j) const;

    // 逻辑索引 → 度量缓存线性偏移（面板内）
    IIndex metric_offset(int panel, int i, int j) const;

    // ---- 坐标生成 ----
    // 预计算全部面板的度量缓存（构造函数中调用；也可显式重算）
    void build_metrics();

    // ---- halo 交换 ----
    // 对定义在单元中心（或交错面）的标量场做跨面板 halo 交换
    //（实现见 halo.cpp；向量场用带旋转的重载）
    void exchange_scalar(std::vector<StateReal>& field) const;

    // 对协变向量分量 (u_cov, v_cov) 做跨面板 halo 交换（含切向旋转）
    void exchange_vector(std::vector<StateReal>& u_cov,
                         std::vector<StateReal>& v_cov) const;

    // ---- 面板连接查询 ----
    // 返回面板 p 在逻辑边 side（0=西,1=东,2=南,3=北）的邻接面板；
    // 返回 -1 表示该边为无邻接（正常情况下立方球 6 面板全部边均有邻接）
    int neighbor_panel(int panel, int side) const;

    // ---- 全球单元统计 ----
    IIndex global_cell_count() const { return 6LL * ncells_ * ncells_; }

private:
    int ncells_ = 48;
    int nhalo_  = 3;
    StateReal radius_ = 6371229.0;

    // 面板内逻辑边长（含 halo）
    int n_ = 0;

    // 度量缓存：每个面板一个扁平数组，长度 n_*n_
    // 布局：metric_cache_[panel * n_*n_ + j * n_ + i]
    std::vector<MetricPoint> metric_cache_;

    // 面板邻接表（[panel][side] -> neighbor panel）
    std::array<std::array<int, 4>, kNumPanels> neighbors_{};

    void init_neighbors();
};

}  // namespace cubed_sph::grid
