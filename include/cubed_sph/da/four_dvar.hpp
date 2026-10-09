// =============================================================================
//  增量 4D-Var 成本函数与极小化（P5）
// =============================================================================
//  增量 4D-Var 的完整代价函数（控制变量形式，δx = B^{1/2} w）：
//      J(w) = ½ wᵀ w
//           + ½ Σ_k (H_k M_k B^{1/2} w - d_k)ᵀ R_k⁻¹ (H_k M_k B^{1/2} w - d_k)
//  其中 d_k = y_k - H_k M_k(x_b) 为创新向量（innovation），M_k 为切线性模型。
//
//  梯度（伴随反推）：
//      ∇_w J = w + B^{1/2ᵀ} Σ_k M_kᵀ H_kᵀ R_k⁻¹ (H_k M_k B^{1/2} w - d_k)
//
//  本模块实现：
//    - FourDVarCostFunction：计算 J(w) 与 ∇J(w)
//    - PreconditionedConjugateGradient：极小化 J（预条件共轭梯度）
//
//  参考文献：
//    Courtier, Thépaut & Hollingsworth (1994), QJRMS 120, 1367–1387
//    Lorenc (2003), QJRMS 129, 3183–3203
//    Rawlins et al. (2007), QJRMS 133, 347–362（混合 4D-Var）
// =============================================================================

#pragma once

#include "cubed_sph/da/observation_operator.hpp"
#include "cubed_sph/da/tlm_ad.hpp"
#include "cubed_sph/io/io.hpp"

#include <vector>
#include <memory>

namespace cubed_sph::dynamics { class DynamicsCore; }
namespace cubed_sph::grid { class VerticalCoordinate; }

namespace cubed_sph::da {

// ---------------------------------------------------------------------------
// 4D-Var 配置
// ---------------------------------------------------------------------------
struct FourDVarConfig {
    int    max_iter = 30;        // 共轭梯度最大迭代
    StateReal grad_tol = 1e-8;   // 梯度收敛容差
    StateReal j_tol = 1e-12;     // 成本函数相对下降容差
};

// ---------------------------------------------------------------------------
// 增量 4D-Var 成本函数
// ---------------------------------------------------------------------------
// 管理一个同化窗内的观测时隙、切线性/伴随模型、背景/观测误差协方差。
// 提供 J(w) 与 ∇J(w) 的计算，供极小化器调用。
// ---------------------------------------------------------------------------
class FourDVarCostFunction {
public:
    // 构造：绑定背景协方差、观测协方差、观测算子、切线性/伴随模型。
    // 观测按时间窗组织为时隙（time slot），每个时隙一组观测。
    FourDVarCostFunction(
        std::shared_ptr<BackgroundCovariance> B,
        std::shared_ptr<ObservationCovariance> R,
        std::shared_ptr<ObservationOperator> H,
        std::shared_ptr<TangentLinearModel> tlm,
        std::shared_ptr<AdjointModel> adm,
        const dynamics::State& x_b,
        StateReal window_len,
        const FourDVarConfig& cfg = {});

    // 计算成本函数值 J(w)
    StateReal evaluate(const dynamics::State& w);

    // 计算梯度 ∇J(w)（写入 grad，原地）
    void gradient(const dynamics::State& w, dynamics::State& grad);

    // 计算创新向量 d_k = y_k - H_k M_k(x_b)（预处理，供 evaluate/gradient 复用）
    void compute_innovations();

    // 设置创新向量 d（由观测值 y 与背景等价量 H(x_b) 计算后注入）
    // 供 FourDVar::run 使用；d 维度须等于 total_nobs()。
    void set_innovation(const std::vector<StateReal>& d);

    // 观测总数与背景
    size_t total_nobs() const { return total_nobs_; }
    const dynamics::State& background() const { return x_b_; }

    // 配置访问器（供极小化器读取收敛参数）
    int max_iter() const { return cfg_.max_iter; }
    StateReal grad_tol() const { return cfg_.grad_tol; }

    // 由控制变量 w 恢复增量 δx = B^{1/2} w
    void increment(const dynamics::State& w, dynamics::State& dx) const;

private:
    std::shared_ptr<BackgroundCovariance> B_;
    std::shared_ptr<ObservationCovariance> R_;
    std::shared_ptr<ObservationOperator> H_;
    std::shared_ptr<TangentLinearModel> tlm_;
    std::shared_ptr<AdjointModel> adm_;

    dynamics::State x_b_;   // 背景（需可拷贝）
    StateReal window_len_;
    FourDVarConfig cfg_;

    size_t total_nobs_ = 0;
    std::vector<StateReal> innovation_;  // 创新向量 d（全部观测拼接）
};

// ---------------------------------------------------------------------------
// 预条件共轭梯度（极小化 J）
// ---------------------------------------------------------------------------
// 增量 4D-Var 中背景项 Jb = ½ wᵀ w 已提供天然预条件（单位阵），
// 故 PCG 的预条件子取单位阵，退化为标准 CG。
// ---------------------------------------------------------------------------
class PreconditionedConjugateGradient {
public:
    // 从初值 w 出发极小化 J，返回迭代次数与最终成本函数值。
    int minimize(FourDVarCostFunction& J,
                 dynamics::State& w,
                 StateReal& final_cost);
};

// ---------------------------------------------------------------------------
// 4D-Var 驱动器（对外高层入口）
// ---------------------------------------------------------------------------
// 串联：构建成本函数 → PCG 极小化 → 输出分析增量 δx = B^{1/2} w*。
// ---------------------------------------------------------------------------
class FourDVar {
public:
    // 执行一次 4D-Var 分析：返回分析状态 x_a = x_b + B^{1/2} w*
    // 输入：背景 x_b、观测集合（含时间/误差/类型）、同化窗长度。
    // 输出：x_a 写入参数 out，返回成本函数迭代信息。
    int run(const dynamics::State& x_b,
            const std::vector<io::Observation>& obs,
            StateReal window_len,
            dynamics::State& x_a,
            StateReal& final_cost);

    // 配置（可调）
    FourDVarConfig cfg;
    StateReal sigma_rho = 0.05;   // 背景误差标准差（密度）
    StateReal sigma_mom = 1.0;
    StateReal sigma_theta = 0.5;

private:
    dynamics::DynamicsCore* core_ = nullptr;          // 外部注入（若需线性化）
    const grid::VerticalCoordinate* vert_ = nullptr;
    StateReal dt_ = 300.0;                            // TL/AD 时间步
};

}  // namespace cubed_sph::da
