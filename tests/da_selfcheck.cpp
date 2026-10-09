// =============================================================================
//  da 模块逻辑自检（不依赖 Config/yaml-cpp，纯核心数学验证）
// =============================================================================
//  验证：
//   1. 观测算子 H 的 TL/AD 严格互为转置（点积一致性）
//   2. 对角背景/观测误差协方差的 B、B^{-1}、B^{1/2}、B^{-1/2} 语义正确
//   3. 增量 4D-Var 成本函数 J(w) 的解析梯度 vs 有限差分梯度
//
//  编译运行：
//   cl /std:c++17 /I include src/da/observation_operator.cpp src/da/tlm_ad.cpp
//       src/da/four_dvar.cpp src/dynamics/state.cpp src/grid/*.cpp 本文件
// =============================================================================

#include "cubed_sph/da/observation_operator.hpp"
#include "cubed_sph/da/tlm_ad.hpp"
#include "cubed_sph/da/four_dvar.hpp"
#include "cubed_sph/io/io.hpp"

#include <cstdio>
#include <cmath>
#include <vector>

using namespace cubed_sph;

int main() {
    // 构造小网格（ncells=6, nhalo=2, nlev=4）
    grid::GridParams gp;
    gp.ncells = 6;
    gp.nhalo = 2;
    grid::CubedSphereGrid g(gp);

    const int nlev = 4;
    dynamics::State x(g, nlev);

    // 填充背景状态（非零，用于 H(x_b) 计算）
    for (IIndex i = 0; i < x.ncell(); ++i) {
        x.rho[i] = 1.2;
        x.rho_theta[i] = 300.0;
        x.theta[i] = 300.0;  // 位温
    }

    // ---- 测试 1：观测算子 TL/AD 转置 ----
    // 构造观测（12 条 TEMP 观测）
    std::vector<io::Observation> obs;
    for (int i = 0; i < 12; ++i) {
        io::Observation o;
        o.value = 295.0 + i * 0.5;
        o.error = 1.0;
        o.lon = 0.1 * i;
        o.lat = 0.2 * i;
        o.type = static_cast<int>(io::ObsType::TEMP);
        obs.push_back(o);
    }

    da::SimpleThetaObservationOperator H(obs, g, nlev);

    // 随机扰动 dx（theta 字段）与 dy（观测空间）
    dynamics::State dx(g, nlev);
    dynamics::State dy_state(g, nlev);
    std::vector<StateReal> dy(obs.size());
    // 简单确定性扰动
    for (IIndex i = 0; i < dx.ncell(); ++i) {
        dx.theta[i] = (i % 7) * 0.01;
    }
    for (size_t i = 0; i < dy.size(); ++i) {
        dy[i] = (i % 5) * 0.1;
    }

    // <H dx, dy>
    std::vector<StateReal> Hdx(obs.size());
    H.tangent_linear(dx, Hdx);
    StateReal lhs = 0.0;
    for (size_t i = 0; i < dy.size(); ++i) lhs += Hdx[i] * dy[i];

    // <dx, H^T dy>
    dynamics::State Htdy(g, nlev);
    std::fill(Htdy.theta.begin(), Htdy.theta.end(), 0.0);
    H.adjoint(dy, Htdy);
    StateReal rhs = 0.0;
    for (IIndex i = 0; i < dx.ncell(); ++i) rhs += dx.theta[i] * Htdy.theta[i];

    StateReal rel1 = std::abs(lhs - rhs) / std::max(std::abs(lhs), std::abs(rhs));
    std::printf("[1] 观测算子 TL/AD 转置: lhs=%.10f rhs=%.10f rel=%.3e %s\n",
                lhs, rhs, rel1, rel1 < 1e-12 ? "PASS" : "FAIL");

    // ---- 测试 2：对角协方差语义 ----
    da::DiagonalBackgroundCovariance B;
    B.sigma_rho = 0.05; B.sigma_mom = 1.0; B.sigma_theta = 0.5;
    dynamics::State v(g, nlev);
    dynamics::State out(g, nlev);
    for (IIndex i = 0; i < v.ncell(); ++i) {
        v.rho[i] = 1.0; v.rho_u[i] = 1.0; v.rho_theta[i] = 1.0;
    }
    // B^{-1/2}(B^{1/2} v) 应恢复 v（rho 字段）
    dynamics::State B12v(g, nlev);
    dynamics::State recovered(g, nlev);
    B.sqrt_B(v, B12v);
    B.sqrt_Binv(B12v, recovered);
    StateReal maxdiff = 0.0;
    for (IIndex i = 0; i < v.ncell(); ++i) {
        maxdiff = std::max(maxdiff, std::abs(recovered.rho[i] - v.rho[i]));
    }
    std::printf("[2] B^{1/2} 与 B^{-1/2} 互逆 (rho): maxdiff=%.3e %s\n",
                maxdiff, maxdiff < 1e-12 ? "PASS" : "FAIL");

    // ---- 测试 3：4D-Var 成本函数梯度（解析 vs 有限差分）----
    auto R = std::make_shared<da::DiagonalObservationCovariance>(obs);
    auto Hp = std::make_shared<da::SimpleThetaObservationOperator>(obs, g, nlev);
    auto Bp = std::make_shared<da::DiagonalBackgroundCovariance>();
    Bp->sigma_rho = B.sigma_rho; Bp->sigma_mom = B.sigma_mom;
    Bp->sigma_theta = B.sigma_theta;

    // TL/AD 用空指针（恒等 M=I，仅测观测项 + 背景项）
    da::FourDVarCostFunction J(Bp, R, Hp, nullptr, nullptr, x, 3600.0);

    // 创新 d = y - H(x_b)
    std::vector<StateReal> hx(obs.size());
    Hp->forward(x, hx);
    std::vector<StateReal> d(obs.size());
    for (size_t i = 0; i < obs.size(); ++i) d[i] = obs[i].value - hx[i];
    J.set_innovation(d);

    // 控制变量 w（非零）
    dynamics::State w(g, nlev);
    for (IIndex i = 0; i < w.ncell(); ++i) {
        w.rho[i] = (i % 3) * 0.001;
        w.rho_theta[i] = (i % 5) * 0.001;
    }

    // 解析梯度
    dynamics::State grad(g, nlev);
    J.gradient(w, grad);

    // 有限差分梯度（对每个控制变量分量）
    const StateReal eps = 1e-6;
    StateReal num = 0.0, den = 0.0;
    for (IIndex i = 0; i < w.ncell(); ++i) {
        // 扰动 rho_theta[i] 分量
        dynamics::State wp(g, nlev); wp = w;
        dynamics::State wm(g, nlev); wm = w;
        wp.rho_theta[i] += eps;
        wm.rho_theta[i] -= eps;
        StateReal Jp = J.evaluate(wp);
        StateReal Jm = J.evaluate(wm);
        StateReal fd = (Jp - Jm) / (2.0 * eps);
        num += (grad.rho_theta[i] - fd) * (grad.rho_theta[i] - fd);
        den += fd * fd;
    }
    StateReal rel3 = std::sqrt(num) / std::sqrt(den);
    std::printf("[3] 4D-Var 梯度解析 vs 有限差分: rel=%.3e %s\n",
                rel3, rel3 < 1e-5 ? "PASS" : "FAIL");

    bool ok = (rel1 < 1e-12) && (maxdiff < 1e-12) && (rel3 < 1e-5);
    std::printf("\n=== %s ===\n", ok ? "ALL PASS" : "SOME FAILED");
    return ok ? 0 : 1;
}
