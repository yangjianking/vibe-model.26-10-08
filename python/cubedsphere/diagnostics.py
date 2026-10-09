"""守恒诊断（纯 Python，与 C++ diagnostics 对应）。

提供全球积分与守恒量计算的参考实现，用于后处理阶段核对 C++ 输出。
采用 Python 的 ``math.fsum``（精确浮点求和，等价于补偿求和）避免大网格
求和误差。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .grid import compute_metric, NUM_PANELS


@dataclass
class ConservationSnapshot:
    """全球守恒量快照。"""

    total_mass: float = 0.0
    total_theta: float = 0.0
    total_energy: float = 0.0
    total_angular_momentum: float = 0.0
    kinetic_energy: float = 0.0
    internal_energy: float = 0.0
    potential_energy: float = 0.0


def global_integral(field, ncells, nlev, z_half, radius=6371229.0) -> float:
    """全球积分 ∫ field dV。

    参数：
        field   : 扁平数组，长度 6×ncells×ncells×nlev（面板, k, j, i 排序）
        ncells  : 每面板单方向单元数
        nlev    : 垂直层数
        z_half  : 半层高度列表（长度 nlev+1），用于层厚
        radius  : 地球半径 [m]
    """
    dxi = 2.0 / ncells
    r2 = radius * radius
    terms = []
    for panel in range(NUM_PANELS):
        for k in range(nlev):
            dz = z_half[k + 1] - z_half[k]
            for j in range(ncells):
                eta = -1.0 + (j + 0.5) * dxi
                for i in range(ncells):
                    xi = -1.0 + (i + 0.5) * dxi
                    m = compute_metric(xi, eta)
                    idx = (((panel * nlev + k) * ncells + j) * ncells + i)
                    terms.append(field[idx] * m.sqrt_g * (r2 * dxi * dxi) * dz)
    return math.fsum(terms)


def relative_drift(initial: ConservationSnapshot,
                   current: ConservationSnapshot) -> dict[str, float]:
    """各守恒量相对初值的相对漂移。"""
    def rel(c, i):
        if abs(i) < 1e-30:
            return 0.0
        return abs(c - i) / abs(i)
    return {
        "mass": rel(current.total_mass, initial.total_mass),
        "theta": rel(current.total_theta, initial.total_theta),
        "energy": rel(current.total_energy, initial.total_energy),
        "angular_momentum": rel(current.total_angular_momentum,
                                initial.total_angular_momentum),
    }
