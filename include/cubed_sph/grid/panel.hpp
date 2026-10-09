// =============================================================================
//  立方球面板枚举（Cubed-Sphere Panels）
// =============================================================================
//  立方球将球面划分为 6 个面板（face），每个面板是一个局部的二维笛卡尔
//  网格。面板编号约定（与 GFDL FV³ / Harris et al. 2021 一致）：
//
//      面板 1–4：赤道带（经度方向依次排列）
//      面板 5：北极帽
//      面板 6：南极帽
//
//  面板间通过共享边（edge）与角点（corner）建立拓扑连接，向量跨面板
//  传递时需做局部坐标系的旋转（见 panel_rotation.hpp）。
//
//  参考文献：
//    Ronchi et al. (1996), J. Comput. Phys. 124, 93–114（"cubed sphere"）
//    Putman & Lin (2007), J. Comput. Phys. 227, 55–78（多种立方球网格）
//    Harris et al. (2021), GFDL Tech. Memo. GFDL2021001 §3（FV³ 立方球）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"

#include <array>

namespace cubed_sph::grid {

// 面板数量
inline constexpr int kNumPanels = 6;

// ---------------------------------------------------------------------------
// 面板 ID（0–5）
// ---------------------------------------------------------------------------
enum class PanelId : int {
    Face1 = 0,   // 赤道带 1（经度 ~ -90° ~ -45° 附近，视具体朝向而定）
    Face2 = 1,   // 赤道带 2
    Face3 = 2,   // 赤道带 3
    Face4 = 3,   // 赤道带 4
    North = 4,   // 北极帽
    South = 5,   // 南极帽
};

inline constexpr PanelId kEquatorialPanels[4] = {
    PanelId::Face1, PanelId::Face2, PanelId::Face3, PanelId::Face4};

inline constexpr PanelId kPolarPanels[2] = {PanelId::North, PanelId::South};

// ---------------------------------------------------------------------------
// 面板边（edge）：立方球有 12 条共享边，每条边连接两个面板
// ---------------------------------------------------------------------------
enum class EdgeId : int {
    // 赤道带面板之间（4 条）
    Face1_Face2 = 0,
    Face2_Face3 = 1,
    Face3_Face4 = 2,
    Face4_Face1 = 3,
    // 赤道带与北极帽（4 条）
    Face1_North = 4,
    Face2_North = 5,
    Face3_North = 6,
    Face4_North = 7,
    // 赤道带与南极帽（4 条）
    Face1_South = 8,
    Face2_South = 9,
    Face3_South = 10,
    Face4_South = 11,
};

inline constexpr int kNumEdges = 12;

// ---------------------------------------------------------------------------
// 面板边方向（用于 halo 交换与向量旋转的符号约定）
// ---------------------------------------------------------------------------
enum class EdgeOrientation : int {
    Same      = 0,   // 邻接面板局部坐标方向一致
    Reversed  = 1,   // 邻接面板局部坐标方向相反
};

// 查询某条边连接的两个面板及方向（在 panel_topology.cpp 中实现）
struct EdgeConnection {
    PanelId a;
    PanelId b;
    EdgeOrientation orientation;
};

// 返回边 e 的连接信息
EdgeConnection edge_connection(EdgeId e);

}  // namespace cubed_sph::grid
