"""立方球网格几何与度量（gnomonic equiangular 投影）。

与 C++ 端 ``include/cubed_sph/grid/metric.hpp`` 对应，纯 Python 实现，
供离线验证与后处理坐标生成。核心公式（Ronchi et al. 1996；Sadourny 1972）：

    g11 = (1 + η²) / r⁴,  g22 = (1 + ξ²) / r⁴,  g12 = -ξ η / r⁴
    √G  = 1 / r³,        r²  = 1 + ξ² + η²

等角坐标 ξ = tan(θ_ξ) ∈ [-1, 1]，θ_ξ ∈ [-π/4, π/4]。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# 面板标识（与 C++ 端 PanelId 对应）
PANELS = ("Face1", "Face2", "Face3", "Face4", "North", "South")
NUM_PANELS = 6

# 地球半径 [m]（与 GridParams 默认一致）
EARTH_RADIUS = 6371229.0


@dataclass(frozen=True)
class MetricPoint:
    """单个网格点的 gnomonic 度量（单位球）。"""

    xi: float
    eta: float
    g11: float
    g12: float
    g22: float
    sqrt_g: float  # √G = 1/r³（单位球）

    @property
    def inv_g11(self) -> float:
        det = self.g11 * self.g22 - self.g12 * self.g12
        return self.g22 / det

    @property
    def inv_g12(self) -> float:
        det = self.g11 * self.g22 - self.g12 * self.g12
        return -self.g12 / det

    @property
    def inv_g22(self) -> float:
        det = self.g11 * self.g22 - self.g12 * self.g12
        return self.g11 / det


def compute_metric(xi: float, eta: float) -> MetricPoint:
    """计算单位球上 (ξ, η) 点的 gnomonic 度量。

    与 C++ ``compute_metric`` 严格对应（不含面板朝向旋转）。
    """
    r2 = 1.0 + xi * xi + eta * eta
    r4 = r2 * r2
    g11 = (1.0 + eta * eta) / r4
    g22 = (1.0 + xi * xi) / r4
    g12 = -xi * eta / r4
    sqrt_g = 1.0 / (r2 * math.sqrt(r2))
    return MetricPoint(xi=xi, eta=eta, g11=g11, g12=g12, g22=g22, sqrt_g=sqrt_g)


# ---------------------------------------------------------------------------
# 面板朝向：将面板局部坐标 (ξ, η) 映射到地理经纬度。
# 6 个面板分别朝向 ±x, ±y, ±z。下面给出标准约定（与 C++ panel_rotation 对应）。
# ---------------------------------------------------------------------------
def panel_lonlat(panel: int, xi: float, eta: float) -> tuple[float, float]:
    """面板局部 (ξ, η) → 经纬度 (lon, lat) [rad]。

    立方球 6 面板分别以 +x, -x, +y, -y, +z, -z 为法向。面板中心
    (ξ=0, η=0) 为法向与单位球的交点。局部切向坐标 (ξ, η) 对应另外两个轴。

    映射约定（法向轴 n，切向轴 t1, t2）：
        Face1(+x):  n=+x, t1=+y, t2=+z
        Face2(-x):  n=-x, t1=-y, t2=+z
        Face3(+y):  n=+y, t1=+x, t2=+z
        Face4(-y):  n=-y, t1=-x, t2=+z
        North(+z):  n=+z, t1=+y, t2=-x
        South(-z):  n=-z, t1=+y, t2=+x
    归一化因子 r = sqrt(1 + ξ² + η²)，切向分量 (ξ, η) 对应 (t1, t2)。
    """
    r = math.sqrt(1.0 + xi * xi + eta * eta)
    # 法向单位分量 = 1/r，切向分量 = ξ/r（沿 t1）、η/r（沿 t2）
    inv_r = 1.0 / r
    xi_n = xi * inv_r
    eta_n = eta * inv_r

    if panel == 0:      # Face1: n=+x, t1=+y, t2=+z
        gx, gy, gz = inv_r, xi_n, eta_n
    elif panel == 1:    # Face2: n=-x, t1=-y, t2=+z
        gx, gy, gz = -inv_r, -xi_n, eta_n
    elif panel == 2:    # Face3: n=+y, t1=+x, t2=+z
        gx, gy, gz = xi_n, inv_r, eta_n
    elif panel == 3:    # Face4: n=-y, t1=-x, t2=+z
        gx, gy, gz = -xi_n, -inv_r, eta_n
    elif panel == 4:    # North: n=+z, t1=+y, t2=-x
        gx, gy, gz = -eta_n, xi_n, inv_r
    else:               # South: n=-z, t1=+y, t2=+x
        gx, gy, gz = eta_n, xi_n, -inv_r

    lon = math.atan2(gy, gx)
    lat = math.asin(max(-1.0, min(1.0, gz)))
    return lon, lat


def panel_coordinates(ncells: int, panel: int) -> tuple[list[float], list[float]]:
    """生成某面板的经纬度坐标数组（单元中心）。

    返回 (lon_flat, lat_flat)，长度 ncells×ncells，按 C++ 的 (j, i) 排序。
    """
    dxi = 2.0 / ncells
    lons: list[float] = []
    lats: list[float] = []
    for j in range(ncells):
        eta = -1.0 + (j + 0.5) * dxi
        for i in range(ncells):
            xi = -1.0 + (i + 0.5) * dxi
            lon, lat = panel_lonlat(panel, xi, eta)
            lons.append(lon)
            lats.append(lat)
    return lons, lats


def unit_sphere_area(ncells: int) -> float:
    """6 面板单位球面积数值积分（用于验证 ≈ 4π）。"""
    dxi = 2.0 / ncells
    total = 0.0
    for panel in range(NUM_PANELS):
        for j in range(ncells):
            eta = -1.0 + (j + 0.5) * dxi
            for i in range(ncells):
                xi = -1.0 + (i + 0.5) * dxi
                m = compute_metric(xi, eta)
                total += m.sqrt_g * dxi * dxi
    return total
