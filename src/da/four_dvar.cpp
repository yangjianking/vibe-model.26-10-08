// =============================================================================
//  增量 4D-Var 成本函数与极小化实现
// =============================================================================
//  实现增量 4D-Var 的完整代价函数 J(w) 与梯度 ∇J(w)，以及预条件共轭梯度
//  极小化。控制变量形式 δx = B^{1/2} w 使背景项退化为 ½ wᵀ w。
//
//  梯度（伴随反推）：
//      ∇_w J = w + B^{1/2ᵀ} Σ_k M_kᵀ H_kᵀ R_k⁻¹ (H_k M_k B^{1/2} w - d_k)
//
//  本实现为单时隙简化（观测集中在同化窗内一个时隙，TEMP 类型），
//  多时隙完整 4D-Var 扩展见 docs/theory/four_dvar.md。
// =============================================================================

#include "cubed_sph/da/four_dvar.hpp"

#include <cmath>
#include <algorithm>

namespace cubed_sph::da {

// =============================================================================
// FourDVarCostFunction
// =============================================================================
FourDVarCostFunction::FourDVarCostFunction(
    std::shared_ptr<BackgroundCovariance> B,
    std::shared_ptr<ObservationCovariance> R,
    std::shared_ptr<ObservationOperator> H,
    std::shared_ptr<TangentLinearModel> tlm,
    std::shared_ptr<AdjointModel> adm,
    const dynamics::State& x_b,
    StateReal window_len,
    const FourDVarConfig& cfg)
    : B_(std::move(B)), R_(std::move(R)), H_(std::move(H)),
      tlm_(std::move(tlm)), adm_(std::move(adm)),
      x_b_(x_b), window_len_(window_len), cfg_(cfg),
      total_nobs_(H_ ? H_->nobs() : 0) {
    innovation_.assign(total_nobs_, 0.0);
}

void FourDVarCostFunction::compute_innovations() {
    // d = y - H(x_b)。完整实现需观测值 y 与背景预报 M(x_b)。
    // 本骨架中创新由 set_innovation() 显式注入（见 FourDVar::run）。
}

void FourDVarCostFunction::set_innovation(const std::vector<StateReal>& d) {
    if (d.size() != total_nobs_) {
        innovation_.assign(total_nobs_, 0.0);
        return;
    }
    innovation_ = d;
}

void FourDVarCostFunction::increment(const dynamics::State& w,
                                     dynamics::State& dx) const {
    // δx = B^{1/2} w
    B_->sqrt_B(w, dx);
}

StateReal FourDVarCostFunction::evaluate(const dynamics::State& w) {
    // J(w) = ½ wᵀ w + ½ (H B^{1/2} w - d)ᵀ R⁻¹ (H B^{1/2} w - d)
    // 背景项 ½ wᵀ w
    StateReal J = 0.5 * state_inner_product(w, w);

    // 增量 δx = B^{1/2} w
    dynamics::State dx(w.grid(), w.nlev());
    increment(w, dx);

    // 观测项：H δx（切线性观测算子），残差 r = H δx - d
    std::vector<StateReal> Hdx(total_nobs_, 0.0);
    H_->tangent_linear(dx, Hdx);
    std::vector<StateReal> r(total_nobs_, 0.0);
    for (size_t i = 0; i < total_nobs_; ++i) {
        r[i] = Hdx[i] - innovation_[i];
    }
    // R⁻¹ r
    std::vector<StateReal> Rinvr(total_nobs_, 0.0);
    R_->apply_Rinv(r, Rinvr);
    // ½ rᵀ R⁻¹ r
    StateReal Jo = 0.0;
    for (size_t i = 0; i < total_nobs_; ++i) {
        Jo += r[i] * Rinvr[i];
    }
    J += 0.5 * Jo;
    return J;
}

void FourDVarCostFunction::gradient(const dynamics::State& w,
                                    dynamics::State& grad) {
    // ∇_w J = w + B^{1/2ᵀ} Hᵀ R⁻¹ (H B^{1/2} w - d)
    // 背景项梯度 = w
    grad = w;

    // 增量 δx = B^{1/2} w
    dynamics::State dx(w.grid(), w.nlev());
    increment(w, dx);

    // H δx
    std::vector<StateReal> Hdx(total_nobs_, 0.0);
    H_->tangent_linear(dx, Hdx);

    // 残差 r = H δx - d，加权 wgt = R⁻¹ r
    std::vector<StateReal> r(total_nobs_, 0.0);
    for (size_t i = 0; i < total_nobs_; ++i) {
        r[i] = Hdx[i] - innovation_[i];
    }
    std::vector<StateReal> wgt(total_nobs_, 0.0);
    R_->apply_Rinv(r, wgt);

    // 观测算子伴随：Hᵀ wgt
    dynamics::State Htw(dx.grid(), dx.nlev());
    // 清零（构造后字段为 0）
    std::fill(Htw.theta.begin(), Htw.theta.end(), 0.0);
    H_->adjoint(wgt, Htw);

    // B^{1/2ᵀ} 作用于 Hᵀ wgt（对角 B 下 B^{1/2} 对称，等于 B^{1/2}）
    dynamics::State B12Htw(dx.grid(), dx.nlev());
    B_->sqrt_B(Htw, B12Htw);

    // grad += B^{1/2} Hᵀ R⁻¹ r
    state_axpy(1.0, B12Htw, grad);
}

// =============================================================================
// PreconditionedConjugateGradient
// =============================================================================
int PreconditionedConjugateGradient::minimize(FourDVarCostFunction& J,
                                              dynamics::State& w,
                                              StateReal& final_cost) {
    // 标准 CG（预条件子为 I，因背景项已是 ½ wᵀ w）
    const StateReal tol = J.grad_tol();

    dynamics::State r(w.grid(), w.nlev());
    dynamics::State p(w.grid(), w.nlev());
    dynamics::State Ap(w.grid(), w.nlev());

    // r0 = -∇J(w0)；p0 = r0
    J.gradient(w, r);
    for (IIndex i = 0; i < r.ncell(); ++i) {
        r.rho[i] = -r.rho[i];
        r.rho_u[i] = -r.rho_u[i];
        r.rho_v[i] = -r.rho_v[i];
        r.rho_w[i] = -r.rho_w[i];
        r.rho_theta[i] = -r.rho_theta[i];
        r.exner[i] = -r.exner[i];
        r.theta[i] = -r.theta[i];
        r.geopotential[i] = -r.geopotential[i];
    }
    p = r;

    StateReal rsold = state_inner_product(r, r);
    final_cost = J.evaluate(w);

    for (int it = 0; it < J.max_iter(); ++it) {
        // 简化：用有限差分 Hessian-vector 积：
        //   ∇²J·p ≈ [∇J(w + ε p) - ∇J(w)] / ε
        const StateReal eps = 1e-7;
        dynamics::State wp(w.grid(), w.nlev());
        wp = w;
        state_axpy(eps, p, wp);
        dynamics::State gwp(w.grid(), w.nlev());
        J.gradient(wp, gwp);
        dynamics::State gw(w.grid(), w.nlev());
        J.gradient(w, gw);
        for (IIndex i = 0; i < w.ncell(); ++i) {
            Ap.rho[i]        = (gwp.rho[i]        - gw.rho[i]) / eps;
            Ap.rho_u[i]      = (gwp.rho_u[i]      - gw.rho_u[i]) / eps;
            Ap.rho_v[i]      = (gwp.rho_v[i]      - gw.rho_v[i]) / eps;
            Ap.rho_w[i]      = (gwp.rho_w[i]      - gw.rho_w[i]) / eps;
            Ap.rho_theta[i]  = (gwp.rho_theta[i]  - gw.rho_theta[i]) / eps;
            Ap.exner[i]      = (gwp.exner[i]      - gw.exner[i]) / eps;
            Ap.theta[i]      = (gwp.theta[i]      - gw.theta[i]) / eps;
            Ap.geopotential[i] = (gwp.geopotential[i] - gw.geopotential[i]) / eps;
        }

        const StateReal pAp = state_inner_product(p, Ap);
        if (std::abs(pAp) < 1e-30) break;
        const StateReal alpha = rsold / pAp;

        state_axpy(alpha, p, w);
        state_axpy(-alpha, Ap, r);

        const StateReal rsnew = state_inner_product(r, r);
        final_cost = J.evaluate(w);
        if (std::sqrt(rsnew) < tol) {
            return it + 1;
        }
        const StateReal beta = rsnew / rsold;
        for (IIndex i = 0; i < p.ncell(); ++i) {
            p.rho[i]        = r.rho[i]        + beta * p.rho[i];
            p.rho_u[i]      = r.rho_u[i]      + beta * p.rho_u[i];
            p.rho_v[i]      = r.rho_v[i]      + beta * p.rho_v[i];
            p.rho_w[i]      = r.rho_w[i]      + beta * p.rho_w[i];
            p.rho_theta[i]  = r.rho_theta[i]  + beta * p.rho_theta[i];
            p.exner[i]      = r.exner[i]      + beta * p.exner[i];
            p.theta[i]      = r.theta[i]      + beta * p.theta[i];
            p.geopotential[i] = r.geopotential[i] + beta * p.geopotential[i];
        }
        rsold = rsnew;
    }
    return J.max_iter();
}

// =============================================================================
// FourDVar
// =============================================================================
int FourDVar::run(const dynamics::State& x_b,
                  const std::vector<io::Observation>& obs,
                  StateReal window_len,
                  dynamics::State& x_a,
                  StateReal& final_cost) {
    // 构建误差协方差与观测算子
    auto B = std::make_shared<DiagonalBackgroundCovariance>();
    B->sigma_rho = sigma_rho;
    B->sigma_mom = sigma_mom;
    B->sigma_theta = sigma_theta;

    auto R = std::make_shared<DiagonalObservationCovariance>(obs);

    const int nlev = x_b.nlev();
    auto H = std::make_shared<SimpleThetaObservationOperator>(
        obs, x_b.grid(), nlev);

    // TL/AD 模型：未注入 core_/vert_ 时用恒等（M = I），仅接口联调。
    // 生产路径需注入 DynamicsCore 与 VerticalCoordinate，构造 HandwrittenTLM/ADM。
    std::shared_ptr<TangentLinearModel> tlm;
    std::shared_ptr<AdjointModel> adm;
    (void)tlm; (void)adm;

    FourDVarCostFunction J(B, R, H, tlm, adm, x_b, window_len, cfg);

    // 创新向量 d = y - H(x_b)
    std::vector<StateReal> hx(obs.size(), 0.0);
    H->forward(x_b, hx);
    std::vector<StateReal> d(obs.size(), 0.0);
    for (size_t i = 0; i < obs.size(); ++i) {
        d[i] = obs[i].value - hx[i];
    }
    J.set_innovation(d);

    // 初值 w = 0（分析 = 背景）
    dynamics::State w(x_b.grid(), nlev);
    std::fill(w.rho.begin(), w.rho.end(), 0.0);
    std::fill(w.rho_u.begin(), w.rho_u.end(), 0.0);
    std::fill(w.rho_v.begin(), w.rho_v.end(), 0.0);
    std::fill(w.rho_w.begin(), w.rho_w.end(), 0.0);
    std::fill(w.rho_theta.begin(), w.rho_theta.end(), 0.0);
    std::fill(w.exner.begin(), w.exner.end(), 0.0);
    std::fill(w.theta.begin(), w.theta.end(), 0.0);
    std::fill(w.geopotential.begin(), w.geopotential.end(), 0.0);

    PreconditionedConjugateGradient pcg;
    final_cost = 0.0;
    int iters = pcg.minimize(J, w, final_cost);

    // 分析增量 δx = B^{1/2} w，x_a = x_b + δx
    dynamics::State dx(x_b.grid(), nlev);
    J.increment(w, dx);
    x_a = x_b;
    state_axpy(1.0, dx, x_a);

    return iters;
}

}  // namespace cubed_sph::da
