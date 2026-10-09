// =============================================================================
//  并行分区与 halo 交换（MPI 进程间 + OpenMP 节点内）
// =============================================================================
//  网格分区：空间填充曲线（Hilbert/Morton）或图划分（METIS）。本阶段提供
//  接口与默认的块分区骨架，完整实现（含 MPI 通信、跨面板 halo 交换、向量
//  旋转）在后续阶段填充，接口保持稳定。
//
//  设计（提示词 P1 并行要求）：
//    - 每 MPI 进程持一个或多个子块；
//    - halo 深度可配置（默认 3，覆盖平流+隐式扩散+嵌套边界缓冲）；
//    - 输出分区质量报告（负载均衡度、通信面/体积比）。
//
//  参考文献：
//    Sagan, H. (1994), Space-Filling Curves（Hilbert 曲线）
//    Karypis & Kumar (1998), SIAM J. Sci. Comput. 20, 359–392（METIS）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"
#include "cubed_sph/grid/cubed_sphere_grid.hpp"

#include <string>
#include <vector>

namespace cubed_sph::parallel {

// ---------------------------------------------------------------------------
// 并行环境（MPI 初始化/终结、rank/size 查询）
// ---------------------------------------------------------------------------
class ParallelContext {
public:
    // 初始化（可选 argc/argv；无 MPI 时退化为单进程）
    static void init(int* argc, char*** argv);
    static void finalize();

    static int rank();
    static int size();
    static bool is_mpi_enabled();
};

// ---------------------------------------------------------------------------
// 域分解：将全球网格单元分配到各 MPI 进程
// ---------------------------------------------------------------------------
class DomainDecomposition {
public:
    DomainDecomposition(const grid::CubedSphereGrid& grid, int nranks);

    // 本进程持有的面板/单元范围
    struct Subdomain {
        int panel;
        int i_begin, i_end;   // 含 halo 的逻辑索引范围
        int j_begin, j_end;
    };

    // 返回本进程（rank）持有的子域列表
    std::vector<Subdomain> subdomains(int rank) const;

    // 分区质量报告（负载均衡度、通信面/体积比）——写入字符串
    std::string quality_report() const;

private:
    const grid::CubedSphereGrid& grid_;
    int nranks_;
};

// ---------------------------------------------------------------------------
// halo 交换（跨进程 + 跨面板）
// ---------------------------------------------------------------------------
// 对定义在面板上的场做 halo 交换。标量场直接拷贝；向量场（协变分量）在
// 跨面板边界需额外做切向旋转（panel_rotation.hpp）。
class HaloExchange {
public:
    explicit HaloExchange(const grid::CubedSphereGrid& grid);

    // 标量场 halo 交换
    void exchange_scalar(std::vector<StateReal>& field);

    // 向量场 halo 交换（含跨面板旋转）
    void exchange_vector(std::vector<StateReal>& u_cov,
                         std::vector<StateReal>& v_cov);

private:
    const grid::CubedSphereGrid& grid_;
};

}  // namespace cubed_sph::parallel
