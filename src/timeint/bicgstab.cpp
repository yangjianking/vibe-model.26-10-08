// =============================================================================
//  BiCGStab 求解器实现（预条件双共轭梯度稳定化）
// =============================================================================
//  BiCGStab 算法（van der Vorst 1992），用于非对称线性系统。
//  伪代码：
//    r0 = b - A x0; r̂0 = r0; ρ0 = α = ω = 1; v = p = 0
//    for j = 1..max_iter:
//      ρj = (r̂0, r_{j-1})
//      β = (ρj/ρ_{j-1})(α/ω)
//      pj = r_{j-1} + β (p_{j-1} - ω v_{j-1})
//      p̂ = M^{-1} pj;  v = A p̂
//      α = ρj / (r̂0, v)
//      s = r_{j-1} - α v
//      若 ‖s‖ < tol：x += α p̂; 返回
//      ŝ = M^{-1} s;  t = A ŝ
//      ω = (t, s)/(t, t)
//      x += α p̂ + ω ŝ
//      r = s - ω t
//      若 ‖r‖ < tol 返回
//
//  参考文献：
//    van der Vorst, H. A. (1992), SIAM J. Sci. Stat. Comput. 13, 631–644
// =============================================================================

#include "cubed_sph/timeint/elliptic_solver.hpp"

#include <cmath>

namespace cubed_sph::timeint {

namespace {

inline StateReal dot(const std::vector<StateReal>& a,
                     const std::vector<StateReal>& b) {
    StateReal s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

inline StateReal norm2(const std::vector<StateReal>& a) {
    return std::sqrt(dot(a, a));
}

}  // namespace

int BiCGStabSolver::solve(const ApplyOp& A, const Precond& M,
                          const std::vector<StateReal>& b,
                          std::vector<StateReal>& x,
                          StateReal tol, int max_iter) {
    const size_t n = b.size();
    if (x.size() != n) x.assign(n, 0.0);

    std::vector<StateReal> r(n), r_hat(n), v(n, 0.0), p(n, 0.0),
        p_hat(n), s(n), s_hat(n), t(n), Ax(n);

    A(x, Ax);
    for (size_t i = 0; i < n; ++i) r[i] = b[i] - Ax[i];
    r_hat = r;

    StateReal rho = 1.0, alpha = 1.0, omega = 1.0;
    const StateReal tol_norm = tol * (norm2(b) + 1e-30);

    for (int iter = 1; iter <= max_iter; ++iter) {
        const StateReal rho_new = dot(r_hat, r);
        if (std::abs(rho_new) < 1e-30) {
            return -iter;  // 分解失败
        }
        if (iter == 1) {
            p = r;
        } else {
            const StateReal beta = (rho_new / rho) * (alpha / omega);
            for (size_t i = 0; i < n; ++i) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
        }

        M(p, p_hat);       // p̂ = M^{-1} p
        A(p_hat, v);       // v = A p̂
        alpha = rho_new / (dot(r_hat, v) + 1e-30);

        for (size_t i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];
        if (norm2(s) < tol_norm) {
            for (size_t i = 0; i < n; ++i) x[i] += alpha * p_hat[i];
            return iter;
        }

        M(s, s_hat);
        A(s_hat, t);
        omega = dot(t, s) / (dot(t, t) + 1e-30);

        for (size_t i = 0; i < n; ++i) {
            x[i] += alpha * p_hat[i] + omega * s_hat[i];
            r[i] = s[i] - omega * t[i];
        }
        if (norm2(r) < tol_norm) return iter;

        if (std::abs(omega) < 1e-30) return -iter;

        rho = rho_new;
    }
    return -max_iter;  // 未收敛
}

}  // namespace cubed_sph::timeint
