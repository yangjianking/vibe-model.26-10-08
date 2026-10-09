// =============================================================================
//  并行环境与域分解实现
// =============================================================================

#include "cubed_sph/parallel/parallel.hpp"

#include <sstream>
#include <string>

#if defined(CS_HAVE_MPI)
#include <mpi.h>
#endif

namespace cubed_sph::parallel {

using grid::kNumPanels;

// ---------------------------------------------------------------------------
// ParallelContext
// ---------------------------------------------------------------------------
void ParallelContext::init(int* argc, char*** argv) {
#if defined(CS_HAVE_MPI)
    MPI_Init(argc, argv);
#else
    (void)argc; (void)argv;  // 无 MPI 时参数未使用
#endif
}

void ParallelContext::finalize() {
#if defined(CS_HAVE_MPI)
    MPI_Finalize();
#endif
}

int ParallelContext::rank() {
#if defined(CS_HAVE_MPI)
    int r; MPI_Comm_rank(MPI_COMM_WORLD, &r); return r;
#else
    return 0;
#endif
}

int ParallelContext::size() {
#if defined(CS_HAVE_MPI)
    int s; MPI_Comm_size(MPI_COMM_WORLD, &s); return s;
#else
    return 1;
#endif
}

bool ParallelContext::is_mpi_enabled() {
#if defined(CS_HAVE_MPI)
    return true;
#else
    return false;
#endif
}

// ---------------------------------------------------------------------------
// DomainDecomposition（块分区骨架：按面板顺序静态均分）
// ---------------------------------------------------------------------------
DomainDecomposition::DomainDecomposition(const grid::CubedSphereGrid& grid,
                                         int nranks)
    : grid_(grid), nranks_(nranks < 1 ? 1 : nranks) {}

std::vector<DomainDecomposition::Subdomain>
DomainDecomposition::subdomains(int rank) const {
    // 骨架：将 6 个面板 × (ncells × ncells) 单元按 rank 顺序静态均分。
    // 完整实现用 Hilbert/Morton 空间填充曲线保证局部性与负载均衡。
    std::vector<Subdomain> out;
    const int n = grid_.panel_n();
    const int per_rank = (kNumPanels * grid_.ncells()) / nranks_;
    const int start = rank * per_rank;
    const int end = (rank == nranks_ - 1) ? kNumPanels * grid_.ncells()
                                          : (rank + 1) * per_rank;

    for (int idx = start; idx < end; ++idx) {
        const int panel = idx / grid_.ncells();
        const int row = idx % grid_.ncells();
        Subdomain s;
        s.panel = panel;
        s.i_begin = grid_.nhalo();
        s.i_end = n - grid_.nhalo();
        s.j_begin = row + grid_.nhalo();
        s.j_end = row + grid_.nhalo() + 1;
        out.push_back(s);
    }
    return out;
}

std::string DomainDecomposition::quality_report() const {
    // 骨架：报告分区总数与面板数。完整负载均衡度/通信面体积比见后续阶段。
    std::ostringstream oss;
    oss << "DomainDecomposition: nranks=" << nranks_
        << ", panels=" << kNumPanels
        << ", ncells/panel=" << grid_.ncells() << "\n";
    oss << "  (load-balance & comm/volume ratio: TODO in later stage)\n";
    return oss.str();
}

// ---------------------------------------------------------------------------
// HaloExchange（骨架：单进程下跨面板拷贝；MPI 通信见后续阶段）
// ---------------------------------------------------------------------------
HaloExchange::HaloExchange(const grid::CubedSphereGrid& grid) : grid_(grid) {}

void HaloExchange::exchange_scalar(std::vector<StateReal>& field) {
    // TODO(parallel): 实现跨进程 halo 交换（MPI_Isend/Irecv + 双缓冲）。
    // 单进程下无跨进程通信，跨面板边界由网格度量连续性保证。
    (void)field;
}

void HaloExchange::exchange_vector(std::vector<StateReal>& u_cov,
                                   std::vector<StateReal>& v_cov) {
    // TODO(parallel): 向量场交换 + 跨面板旋转（panel_rotation）。
    (void)u_cov; (void)v_cov;
}

}  // namespace cubed_sph::parallel
