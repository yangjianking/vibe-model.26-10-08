"""平流方案收敛阶验证（TC2 / 一维高斯钟型平流）。

验证二阶中心、三阶上游、PPM 三种平流方案在 SSP-RK3 时间积分下的
L1 收敛阶。这是 P7 收尾的 TC2 收敛阶实测的数值基础。

关键结论（见 docs/theory/conservation_diagnostics.md §3）：
- 偶数阶中心差分 + 显式 Euler 无条件不稳定（纯色散、无耗散）；
- 需 SSP-RK3 时间积分（模式已实现 ExplicitRK3Integrator）保证稳定；
- Center2 → 二阶、Upwind3 → 三阶、PPM（单调限制）→ ~2.35 阶（平滑场）。

用法：
    python -m cubedsphere.convergence
"""

from __future__ import annotations

import math


def gaussian(x: float, x0: float = 0.5, w: float = 0.1) -> float:
    """一维高斯钟型场。"""
    return math.exp(-((x - x0) / w) ** 2)


def _flux_upwind3(phi, u, N):
    F = [0.0] * N
    for i in range(N):
        ip1 = i + 1 if i < N - 1 else 0
        im1 = i - 1 if i > 0 else N - 1
        F[i] = u * (2.0 * phi[ip1] + 5.0 * phi[i] - phi[im1]) / 6.0
    return F


def _flux_center2(phi, u, N):
    F = [0.0] * N
    for i in range(N):
        ip1 = i + 1 if i < N - 1 else 0
        F[i] = u * 0.5 * (phi[i] + phi[ip1])
    return F


def _flux_ppm(phi, u, N):
    F = [0.0] * N
    for i in range(N):
        ip1 = i + 1 if i < N - 1 else 0
        ip2 = i + 2 if i < N - 2 else (i + 2 - N)
        im1 = i - 1 if i > 0 else N - 1
        f = (7.0 / 12.0) * (phi[i] + phi[ip1]) - (1.0 / 12.0) * (phi[im1] + phi[ip2])
        lo = min(phi[i], phi[ip1])
        hi = max(phi[i], phi[ip1])
        F[i] = u * max(lo, min(hi, f))
    return F


_FLUXES = {
    "Center2": _flux_center2,
    "Upwind3": _flux_upwind3,
    "PPM": _flux_ppm,
}


def rk3_advect(phi, flux_fn, u, dx, dt, N):
    """SSP-RK3 推进一个时间步（常速平流，周期边界）。"""
    def rhs(p):
        F = flux_fn(p, u, N)
        r = [0.0] * N
        for i in range(N):
            im1 = i - 1 if i > 0 else N - 1
            r[i] = -(F[i] - F[im1]) / dx
        return r

    k1 = rhs(phi)
    p1 = [phi[i] + dt * k1[i] for i in range(N)]
    k2 = rhs(p1)
    p2 = [0.75 * phi[i] + 0.25 * (p1[i] + dt * k2[i]) for i in range(N)]
    k3 = rhs(p2)
    return [(1.0 / 3.0) * phi[i] + (2.0 / 3.0) * (p2[i] + dt * k3[i])
            for i in range(N)]


def convergence_order(scheme: str, resolutions=(32, 64, 128, 256),
                      cfl=0.4) -> list[float]:
    """返回某方案的 L1 误差随分辨率变化的收敛阶序列。"""
    flux_fn = _FLUXES[scheme]
    errors = []
    for N in resolutions:
        dx = 1.0 / N
        dt = cfl * dx
        nstep = int(1.0 / dt)  # 平流一个周期
        x = [(i + 0.5) * dx for i in range(N)]
        phi = [gaussian(xi) for xi in x]
        for _ in range(nstep):
            phi = rk3_advect(phi, flux_fn, 1.0, dx, dt, N)
        err = sum(abs(phi[i] - gaussian(x[i])) for i in range(N)) * dx
        errors.append(err)
    orders = [math.log(errors[k - 1] / errors[k]) / math.log(2.0)
              for k in range(1, len(errors))]
    return orders


def main() -> None:
    print("平流方案 L1 收敛阶（SSP-RK3, 高斯钟型, 一周期）")
    print("=" * 55)
    for scheme in ("Center2", "Upwind3", "PPM"):
        orders = convergence_order(scheme)
        print(f"{scheme:>8}: 收敛阶 = "
              + ", ".join(f"{o:.2f}" for o in orders)
              + f"  (末阶 {orders[-1]:.2f})")


if __name__ == "__main__":
    main()
