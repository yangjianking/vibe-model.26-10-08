// =============================================================================
//  观测算子与误差协方差（P5 4D-Var 基础）
// =============================================================================
//  增量 4D-Var 的成本函数为：
//      J(δx) = ½ δxᵀ B⁻¹ δx
//            + ½ Σ_k (H_k M_k δx - d_k)ᵀ R_k⁻¹ (H_k M_k δx - d_k)
//  其中：
//      δx   —— 初始时刻的控制增量向量（相对背景 x_b）
//      B    —— 背景误差协方差矩阵
//      M_k  —— 从初始时刻到观测时刻 t_k 的非线性/切线性预报模型
//      H_k  —— 观测算子（模式空间 → 观测空间）
//      d_k  —— 创新向量（innovation）：d_k = y_k - H_k M_k(x_b)
//      R_k  —— 观测误差协方差矩阵
//
//  本模块提供三个正交接口：
//    - ObservationOperator（H）：模式状态 → 观测等价量，及其切线性/伴随
//    - BackgroundCovariance（B）：背景误差协方差的应用 B v、B⁻¹ v
//    - ObservationCovariance（R）：观测误差协方差的应用 R v、R⁻¹ v
//
//  增量形式（Courtier et al. 1994）的核心思想：不用完整的 B，而是通过
//  控制变量变换 δx = B^{1/2} w 隐式施加背景约束，避免显式构造稠密 B。
//
//  参考文献：
//    Courtier, Thépaut & Hollingsworth (1994), QJRMS 120, 1367–1387
//    Lorenc (2003), QJRMS 129, 3183–3203（增量 4D-Var）
//    Bannister (2008), QJRMS 134, 1951–1970（背景误差协方差建模）
// =============================================================================

#pragma once

#include "cubed_sph/dynamics/state.hpp"
#include "cubed_sph/io/io.hpp"

#include <vector>
#include <functional>

namespace cubed_sph::da {

// ---------------------------------------------------------------------------
// 观测算子接口 H：模式状态 → 观测等价量
// ---------------------------------------------------------------------------
// 观测算子将模式状态插值到观测位置（水平插值 + 垂直插值 + 变量变换）。
// 本接口以"观测向量维度"为输入输出粒度，避免直接耦合模式网格。
// ---------------------------------------------------------------------------
class ObservationOperator {
public:
    virtual ~ObservationOperator() = default;

    // 非线性观测算子：给定模式状态，计算观测等价量 H(x)
    // 输入 state：模式状态；输出 hx：观测等价量（维度 = 观测数 n_obs）
    virtual void forward(const dynamics::State& state,
                         std::vector<StateReal>& hx) const = 0;

    // 切线性观测算子：H δx（对模式状态扰动的线性化）
    // 输入 dx：模式状态扰动；输出 hdx：观测空间扰动
    virtual void tangent_linear(const dynamics::State& dx,
                                std::vector<StateReal>& hdx) const = 0;

    // 伴随观测算子：Hᵀ δy（观测空间扰动 → 模式状态扰动）
    // 输入 dy：观测空间扰动；输出 atdx：模式状态扰动（累加）
    virtual void adjoint(const std::vector<StateReal>& dy,
                         dynamics::State& atdx) const = 0;

    // 观测维度（返回当前算子对应的观测数）
    virtual size_t nobs() const = 0;
};

// ---------------------------------------------------------------------------
// 背景误差协方差 B 接口
// ---------------------------------------------------------------------------
// B 为巨型矩阵，绝不显式构造。本接口仅提供矩阵-向量积：
//    - apply_B(v)     = B v      （用于随机扰动生成、分析误差诊断）
//    - apply_Binv(v)  = B⁻¹ v    （用于成本函数背景项）
//    - sqrt_B(w)      = B^{1/2} w（控制变量变换 δx = B^{1/2} w）
// 默认实现采用对角 B（方差场），完整实现需背景误差相关模型（球谐/小波）。
// ---------------------------------------------------------------------------
class BackgroundCovariance {
public:
    virtual ~BackgroundCovariance() = default;

    // B v
    virtual void apply_B(const dynamics::State& v, dynamics::State& Bv) const = 0;
    // B⁻¹ v
    virtual void apply_Binv(const dynamics::State& v,
                            dynamics::State& Binvv) const = 0;
    // B^{1/2} w（控制变量变换正变换）
    virtual void sqrt_B(const dynamics::State& w, dynamics::State& B12w) const = 0;
    // B^{-1/2} v（控制变量变换逆变换）
    virtual void sqrt_Binv(const dynamics::State& v,
                           dynamics::State& Bmin12v) const = 0;
};

// ---------------------------------------------------------------------------
// 对角背景误差协方差（默认实现）
// ---------------------------------------------------------------------------
// 每个网格点（变量）独立方差 σ_b²，B、B⁻¹、B^{1/2}、B^{-1/2} 均退化为逐点缩放。
// 方差场由外部注入（可来自气候态、NMC 方法或静态猜测）。
// ---------------------------------------------------------------------------
class DiagonalBackgroundCovariance : public BackgroundCovariance {
public:
    // 各预后变量的背景误差标准差（按变量分别设置）
    StateReal sigma_rho = 0.05;   // 密度背景误差标准差（相对比例，可调）
    StateReal sigma_mom = 1.0;    // 动量密度背景误差标准差
    StateReal sigma_theta = 0.5;  // 位温密度背景误差标准差

    void apply_B(const dynamics::State& v, dynamics::State& Bv) const override;
    void apply_Binv(const dynamics::State& v, dynamics::State& Binvv) const override;
    void sqrt_B(const dynamics::State& w, dynamics::State& B12w) const override;
    void sqrt_Binv(const dynamics::State& v, dynamics::State& Bmin12v) const override;

private:
    // 每个字段的方差 σ_b²（按字段返回）
    StateReal var_for_field(int field_id) const;
};

// ---------------------------------------------------------------------------
// 观测误差协方差 R 接口
// ---------------------------------------------------------------------------
// R 通常为对角阵（各观测独立误差），也支持块对角（同类观测相关）。
// 默认实现为对角：R v 与 R⁻¹ v 均为逐观测缩放。
// ---------------------------------------------------------------------------
class ObservationCovariance {
public:
    virtual ~ObservationCovariance() = default;

    // R v
    virtual void apply_R(const std::vector<StateReal>& v,
                         std::vector<StateReal>& Rv) const = 0;
    // R⁻¹ v
    virtual void apply_Rinv(const std::vector<StateReal>& v,
                            std::vector<StateReal>& Rinvv) const = 0;

    virtual size_t nobs() const = 0;
};

// ---------------------------------------------------------------------------
// 对角观测误差协方差（默认实现）
// ---------------------------------------------------------------------------
// 每条观测独立，误差标准差来自 Observation::error。
// ---------------------------------------------------------------------------
class DiagonalObservationCovariance : public ObservationCovariance {
public:
    // 由观测集合构造（提取每条观测的误差标准差 σ_o）
    explicit DiagonalObservationCovariance(const std::vector<io::Observation>& obs);

    void apply_R(const std::vector<StateReal>& v,
                 std::vector<StateReal>& Rv) const override;
    void apply_Rinv(const std::vector<StateReal>& v,
                    std::vector<StateReal>& Rinvv) const override;

    size_t nobs() const override { return var_.size(); }

private:
    std::vector<StateReal> var_;  // σ_o²（每条观测）
};

// ---------------------------------------------------------------------------
// 简单水平/垂直插值观测算子（示例实现）
// ---------------------------------------------------------------------------
// 演示观测算子接口：将模式状态按经纬度水平插值（最近邻/双线性占位）并
// 做变量变换到观测等价量。此处以"位温 θ"为例（观测 TEMP → 模式 θ）。
// 完整实现需接入网格度量与垂直坐标，此处给出可编译、可点积检验的骨架。
// ---------------------------------------------------------------------------
class SimpleThetaObservationOperator : public ObservationOperator {
public:
    // 由观测集合构造，绑定模式网格（用于坐标插值）
    SimpleThetaObservationOperator(const std::vector<io::Observation>& obs,
                                   const grid::CubedSphereGrid& grid,
                                   int nlev);

    void forward(const dynamics::State& state,
                 std::vector<StateReal>& hx) const override;
    void tangent_linear(const dynamics::State& dx,
                        std::vector<StateReal>& hdx) const override;
    void adjoint(const std::vector<StateReal>& dy,
                 dynamics::State& atdx) const override;

    size_t nobs() const override { return nobs_; }

private:
    std::vector<io::Observation> obs_;
    const grid::CubedSphereGrid& grid_;
    int nlev_;
    size_t nobs_ = 0;
};

}  // namespace cubed_sph::da
