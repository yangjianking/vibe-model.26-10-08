// =============================================================================
//  Helmholtz / 椭圆求解器接口（P2 第 4 条）
// =============================================================================
//  半隐式时间推进的核心是求解三维椭圆方程：
//    (I - β² Δt² c_s² ∇²) δπ' = RHS
//
//  求解策略：预条件 Krylov 子空间法（BiCGStab/GMRES）+ 几何多重网格预条件。
//  立方球网格上椭圆算子具有各向异性系数（度量张量 g_ij），需特殊处理；
//  垂直方向强耦合用线松弛/ADI；matrix-free 实现，算子以 stencil 回调暴露。
//
//  本阶段定义抽象接口（SolverBase），提供 Jacobi 预条件 + BiCGStab 骨架，
//  多重网格预条件作为可扩展策略接入。
//
//  参考文献：
//    Buckeridge & Scheichl (2010), QJRMS 136, 2723–2735（球面多重网格）
//    Dedner et al. (2016), Geosci. Model Dev. 9, 2223–2241（GungHo 椭圆）
//    Knoll & Keyes (2004), J. Comput. Phys. 193, 357–397（JFNK）
// =============================================================================

#pragma once

#include "cubed_sph/common/types.hpp"
#include "cubed_sph/common/config.hpp"

#include <vector>
#include <functional>

namespace cubed_sph::timeint {

// ---------------------------------------------------------------------------
// 椭圆求解器抽象接口
// ---------------------------------------------------------------------------
class EllipticSolver {
public:
    virtual ~EllipticSolver() = default;

    // 求解 A x = b。矩阵 A 以 matrix-free 算子回调给出（A(apply) 作用于
    // 向量），预条件器 M^{-1} 可选。返回迭代次数。
    using ApplyOp = std::function<void(const std::vector<StateReal>&,
                                       std::vector<StateReal>&)>;
    using Precond = std::function<void(const std::vector<StateReal>&,
                                       std::vector<StateReal>&)>;

    virtual int solve(const ApplyOp& A, const Precond& M,
                      const std::vector<StateReal>& b,
                      std::vector<StateReal>& x,
                      StateReal tol = 1e-10, int max_iter = 100) = 0;

    virtual const char* name() const = 0;
};

// ---------------------------------------------------------------------------
// Jacobi 预条件（对角近似，骨架）
// ---------------------------------------------------------------------------
class JacobiPreconditioner {
public:
    explicit JacobiPreconditioner(const std::vector<StateReal>& diag)
        : diag_(diag) {}

    void operator()(const std::vector<StateReal>& r,
                    std::vector<StateReal>& z) const {
        for (size_t i = 0; i < r.size(); ++i) {
            z[i] = r[i] / (diag_[i] + 1e-30);
        }
    }

private:
    std::vector<StateReal> diag_;
};

// ---------------------------------------------------------------------------
// BiCGStab 求解器（预条件双共轭梯度稳定化）
// ---------------------------------------------------------------------------
// 实现见 src/timeint/bicgstab.cpp。适合非对称椭圆算子（半隐式可压缩方程
// 的 Helmholtz 算子一般非对称，因含平流修正项）。
// ---------------------------------------------------------------------------
class BiCGStabSolver : public EllipticSolver {
public:
    int solve(const ApplyOp& A, const Precond& M,
              const std::vector<StateReal>& b,
              std::vector<StateReal>& x,
              StateReal tol = 1e-10, int max_iter = 100) override;

    const char* name() const override { return "BiCGStab"; }
};

}  // namespace cubed_sph::timeint
