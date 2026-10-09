// =============================================================================
//  立方球网格实现：度量构造 + 邻接拓扑
// =============================================================================
//  核心：gnomonic equiangular 投影的度量张量与面板邻接关系。
//
//  gnomonic 投影（Ronchi et al. 1996；Putman & Lin 2007）：
//   面板局部等角坐标 (ξ, η) ∈ [-1, 1]，ξ = tan(θ_ξ)，θ_ξ ∈ [-π/4, π/4]。
//   映射到球面单位向量
//     x = (1, ξ, η) / r,   r = sqrt(1 + ξ² + η²)
//   度量张量（协变分量）：
//     g11 = (1 + η²) / r⁴
//     g22 = (1 + ξ²) / r⁴
//     g12 = - ξ η / r⁴
//   Jacobian 面积元： √G = sqrt(det g) = 1 / r³
//
// 完整推导（含朝向矩阵与经纬度反算）见
//   docs/theory/cubed_sphere_metric.md
// =============================================================================

#include "cubed_sph/grid/cubed_sphere_grid.hpp"

#include <cmath>
#include <array>

namespace cubed_sph::grid {

namespace {

// 面板朝向（面板中心方向在笛卡尔系下的单位向量）
// 约定（与 Harris et al. 2021 一致的常见朝向）：
//   Face1 = +x, Face2 = -y, Face3 = -x, Face4 = +y, North = +z, South = -z
constexpr std::array<std::array<StateReal, 3>, 6> kPanelCenter = {{
    {{ 1.0, 0.0, 0.0}},   // Face1
    {{ 0.0, -1.0, 0.0}},  // Face2
    {{-1.0, 0.0, 0.0}},   // Face3
    {{ 0.0, 1.0, 0.0}},   // Face4
    {{ 0.0, 0.0, 1.0}},   // North
    {{ 0.0, 0.0, -1.0}},  // South
}};

// 面板局部基向量 e1（ξ 方向）与 e2（η 方向），基于面板朝向构造。
// 这里为每个面板给出其切平面内的两个正交单位向量（右手系）。
constexpr std::array<std::array<StateReal, 3>, 6> kPanelE1 = {{
    {{ 0.0, 1.0, 0.0}},   // Face1: e1 = +y
    {{-1.0, 0.0, 0.0}},   // Face2: e1 = -x
    {{ 0.0, -1.0, 0.0}},  // Face3: e1 = -y
    {{ 1.0, 0.0, 0.0}},   // Face4: e1 = +x
    {{ 1.0, 0.0, 0.0}},   // North: e1 = +x
    {{ 1.0, 0.0, 0.0}},   // South: e1 = +x
}};

constexpr std::array<std::array<StateReal, 3>, 6> kPanelE2 = {{
    {{ 0.0, 0.0, 1.0}},   // Face1: e2 = +z
    {{ 0.0, 0.0, 1.0}},   // Face2: e2 = +z
    {{ 0.0, 0.0, 1.0}},   // Face3: e2 = +z
    {{ 0.0, 0.0, 1.0}},   // Face4: e2 = +z
    {{ 0.0, 1.0, 0.0}},   // North: e2 = +y
    {{ 0.0, -1.0, 0.0}},  // South: e2 = -y
}};

// 由面板局部坐标 (ξ, η) 反算球面笛卡尔单位向量
inline void gnomonic_to_cart(int panel, StateReal xi, StateReal eta,
                             StateReal& x, StateReal& y, StateReal& z) {
    const StateReal r2 = 1.0 + xi * xi + eta * eta;
    const StateReal r = std::sqrt(r2);
    const StateReal c = kPanelCenter[panel][0];
    const StateReal c1 = kPanelE1[panel][0];
    const StateReal c2 = kPanelE2[panel][0];
    const StateReal e1x = kPanelE1[panel][0], e1y = kPanelE1[panel][1], e1z = kPanelE1[panel][2];
    const StateReal e2x = kPanelE2[panel][0], e2y = kPanelE2[panel][1], e2z = kPanelE2[panel][2];
    (void)c; (void)c1; (void)c2;
    // 球面点 = (n̂ + ξ·e1 + η·e2) / r
    x = (kPanelCenter[panel][0] + xi * e1x + eta * e2x) / r;
    y = (kPanelCenter[panel][1] + xi * e1y + eta * e2y) / r;
    z = (kPanelCenter[panel][2] + xi * e1z + eta * e2z) / r;
}

}  // namespace

// ---------------------------------------------------------------------------
// 度量计算（gnomonic equiangular）
// ---------------------------------------------------------------------------
MetricPoint compute_metric(StateReal xi, StateReal eta, int panel) {
    MetricPoint m;
    m.xi = xi;
    m.eta = eta;

    const StateReal r2 = 1.0 + xi * xi + eta * eta;
    const StateReal r  = std::sqrt(r2);
    const StateReal r4 = r2 * r2;

    // 协变度量张量（Putman & Lin 2007, Eq. A1-A6）
    m.g11 = (1.0 + eta * eta) / r4;
    m.g22 = (1.0 + xi * xi) / r4;
    m.g12 = -xi * eta / r4;

    // 行列式 det(g) = 1/r^6，因此 √G = 1/r³
    m.sqrtG = 1.0 / (r2 * r);

    // 逆度量（反变分量）：对 2x2 对称矩阵求逆
    const StateReal det = m.g11 * m.g22 - m.g12 * m.g12;
    m.inv_g11 =  m.g22 / det;
    m.inv_g22 =  m.g11 / det;
    m.inv_g12 = -m.g12 / det;

    // 球面笛卡尔坐标与经纬度
    gnomonic_to_cart(panel, xi, eta, m.cart_x, m.cart_y, m.cart_z);
    m.normal_x = m.cart_x;
    m.normal_y = m.cart_y;
    m.normal_z = m.cart_z;

    // 经纬度（弧度）
    m.lat = std::asin(m.cart_z);
    m.lon = std::atan2(m.cart_y, m.cart_x);

    return m;
}

// ---------------------------------------------------------------------------
// 邻接表初始化
// ---------------------------------------------------------------------------
void CubedSphereGrid::init_neighbors() {
    // 面板间连接（约定见 panel.hpp 的 EdgeId）
    // side 顺序：0=西(-ξ), 1=东(+ξ), 2=南(-η), 3=北(+η)
    // 以下邻接基于 kPanelCenter/kPanelE1/kPanelE2 的朝向约定手工标注，
    // 完整旋转矩阵在 panel_rotation.cpp 中离线构造。
    //
    // 赤道带首尾相接：Face1-Face2-Face3-Face4-Face1
    neighbors_[0][1] = 1;  // Face1 东 -> Face2
    neighbors_[1][0] = 0;  // Face2 西 -> Face1
    neighbors_[1][1] = 2;  // Face2 东 -> Face3
    neighbors_[2][0] = 1;  // Face3 西 -> Face2
    neighbors_[2][1] = 3;  // Face3 东 -> Face4
    neighbors_[3][0] = 2;  // Face4 西 -> Face3
    neighbors_[3][1] = 0;  // Face4 东 -> Face1
    neighbors_[0][0] = 3;  // Face1 西 -> Face4

    // 赤道带北边 -> 北极帽（北极帽四边依次接 Face1/2/3/4）
    neighbors_[0][3] = 4;  // Face1 北 -> North
    neighbors_[1][3] = 4;  // Face2 北 -> North
    neighbors_[2][3] = 4;  // Face3 北 -> North
    neighbors_[3][3] = 4;  // Face4 北 -> North
    // 北极帽四条边（约定北帽西/东/南/北依次接 Face1/2/3/4）
    neighbors_[4][0] = 1;  // North 西 -> Face2（近似，完整见拓扑文档）
    neighbors_[4][1] = 3;  // North 东 -> Face4
    neighbors_[4][2] = 2;  // North 南 -> Face3
    neighbors_[4][3] = 0;  // North 北 -> Face1

    // 赤道带南边 -> 南极帽
    neighbors_[0][2] = 5;  // Face1 南 -> South
    neighbors_[1][2] = 5;  // Face2 南 -> South
    neighbors_[2][2] = 5;  // Face3 南 -> South
    neighbors_[3][2] = 5;  // Face4 南 -> South
    neighbors_[5][0] = 1;  // South 西 -> Face2
    neighbors_[5][1] = 3;  // South 东 -> Face4
    neighbors_[5][2] = 2;  // South 南 -> Face3
    neighbors_[5][3] = 0;  // South 北 -> Face1

    // 初始化为 -1 表示无邻接的占位（此处 6 面板全部边均有邻接）
    for (auto& row : neighbors_) {
        for (auto& v : row) {
            if (v == 0 && &v == &row[0]) continue;  // 占位保护
            if (v < 0) v = -1;
        }
    }
    // 修正：把未显式赋值的 -1 保留（理论上不会出现）
}

// ---------------------------------------------------------------------------
// 构造函数
// ---------------------------------------------------------------------------
CubedSphereGrid::CubedSphereGrid(const GridParams& params)
    : ncells_(params.ncells),
      nhalo_(params.nhalo),
      radius_(params.radius),
      n_(params.ncells + 2 * params.nhalo) {
    for (auto& row : neighbors_) row.fill(-1);
    init_neighbors();
    metric_cache_.resize(static_cast<size_t>(kNumPanels) * n_ * n_);
    build_metrics();
}

// ---------------------------------------------------------------------------
// 度量构造
// ---------------------------------------------------------------------------
void CubedSphereGrid::build_metrics() {
    // 等角坐标范围 [-1, 1]（ξ = tan θ, θ ∈ [-π/4, π/4]），单元中心坐标：
    //   ξ_i = -1 + (i - nhalo + 0.5) * dxi   （i 为含 halo 的逻辑索引）
    // 其中 dxi = 2 / ncells。
    const StateReal dxi = 2.0 / static_cast<StateReal>(ncells_);
    for (int p = 0; p < kNumPanels; ++p) {
        for (int j = 0; j < n_; ++j) {
            for (int i = 0; i < n_; ++i) {
                const StateReal xi = -1.0 +
                    (static_cast<StateReal>(i - nhalo_) + 0.5) * dxi;
                const StateReal eta = -1.0 +
                    (static_cast<StateReal>(j - nhalo_) + 0.5) * dxi;
                metric_cache_[metric_offset(p, i, j)] = compute_metric(xi, eta, p);
            }
        }
    }
}

IIndex CubedSphereGrid::metric_offset(int panel, int i, int j) const {
    return static_cast<IIndex>(panel) * n_ * n_ + static_cast<IIndex>(j) * n_ + i;
}

const MetricPoint& CubedSphereGrid::metric(int panel, int i, int j) const {
    return metric_cache_[metric_offset(panel, i, j)];
}

int CubedSphereGrid::neighbor_panel(int panel, int side) const {
    if (panel < 0 || panel >= kNumPanels || side < 0 || side >= 4) return -1;
    return neighbors_[panel][side];
}

// ---------------------------------------------------------------------------
// 标量 halo 交换（占位：单进程下无需跨进程通信，跨面板边界填充在
// halo.cpp 中实现；本阶段提供接口与串行跨面板拷贝骨架）
// ---------------------------------------------------------------------------
void CubedSphereGrid::exchange_scalar(std::vector<StateReal>& field) const {
    // 单进程串行场景：字段布局与 metric_cache_ 一致，跨面板边界通过
    // 邻接表拷贝。多进程 MPI 版见 src/parallel/halo_exchange.cpp。
    // 本阶段仅做数据布局校验占位，完整实现留给 P4（并行 I/O）联动。
    (void)field;
}

void CubedSphereGrid::exchange_vector(std::vector<StateReal>& u_cov,
                                      std::vector<StateReal>& v_cov) const {
    // 向量场交换需额外做切向旋转（panel_rotation），完整实现见 halo.cpp。
    (void)u_cov;
    (void)v_cov;
}

}  // namespace cubed_sph::grid
