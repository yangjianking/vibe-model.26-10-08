"""标准算例解析解（Williamson TC2 等）。

与 C++ ``src/testcases/test_cases.cpp`` 对应，用于后处理阶段的误差计算与
收敛阶验证。纯 Python，无重依赖。
"""

from __future__ import annotations

import math

EARTH_RADIUS = 6371229.0


def _great_circle_dist(lon1, lat1, lon2, lat2) -> float:
    """两点间大圆角距离 [rad]。"""
    cosd = (math.sin(lat1) * math.sin(lat2)
            + math.cos(lat1) * math.cos(lat2) * math.cos(lon1 - lon2))
    return math.acos(max(-1.0, min(1.0, cosd)))


class WilliamsonTC2:
    """稳态定常平流（余弦钟型场被定常风场平流）。

    解析解为钟型场绕地球刚体平移，12 天回到原点。
    """

    def __init__(self, alpha=math.pi / 4.0,
                 lambda_c=-3.0 * math.pi / 4.0,
                 phi_c=0.0,
                 bell_radius=math.pi / 3.0):
        self.alpha = alpha
        self.lambda_c = lambda_c
        self.phi_c = phi_c
        self.bell_radius = bell_radius
        # 平流角速度 [rad/s]
        self.u0 = 2.0 * math.pi * EARTH_RADIUS / (12.0 * 86400.0)
        self.omega_adv = self.u0 / EARTH_RADIUS

    def analytic(self, lon: float, lat: float, t: float) -> float:
        """t 秒后的钟型场值（0–1 归一化）。"""
        omega_t = self.omega_adv * t
        # 纯纬向平流分量（α 的经向投影）；骨架实现与 C++ 一致
        lon_prime = lon - omega_t * math.sin(self.alpha)
        lat_prime = lat
        d = _great_circle_dist(lon_prime, lat_prime, self.lambda_c, self.phi_c)
        if d < self.bell_radius:
            return 0.5 * (1.0 + math.cos(math.pi * d / self.bell_radius))
        return 0.0


class HeldSuarez:
    """Held–Suarez (1994) 气候态强迫的平衡温度与摩擦。"""

    KAPPA = 287.05 / 1004.5  # R_d / c_p
    Z_TOP = 30000.0
    SIGMA_B = 0.7

    @classmethod
    def equilibrium_temperature(cls, z: float, lat: float) -> float:
        """平衡温度剖面 T_eq(φ, σ) [K]，σ = p/p_s ≈ 1 - z/z_top。"""
        sigma = 1.0 - z / cls.Z_TOP
        s2 = math.sin(lat) ** 2
        c2 = math.cos(lat) ** 2
        tropo = 0.0
        if sigma > cls.SIGMA_B:
            tropo = (sigma - cls.SIGMA_B) / (1.0 - cls.SIGMA_B)
        t = (315.0 - 60.0 * s2 - 10.0 * c2 * tropo) * sigma ** cls.KAPPA
        return max(200.0, t)

    @classmethod
    def rayleigh_friction(cls, z: float) -> float:
        """边界层瑞利摩擦系数 [s⁻¹]，σ > σ_b 时非零。"""
        sigma = 1.0 - z / cls.Z_TOP
        if sigma > cls.SIGMA_B:
            return (1.0 / 86400.0) * (sigma - cls.SIGMA_B) / (1.0 - cls.SIGMA_B)
        return 0.0
