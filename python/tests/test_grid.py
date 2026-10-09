"""核心几何与度量单元测试。"""

import math
import sys
import os

# 使 tests 可导入 cubedsphere 包（无需安装）
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from cubedsphere.grid import (
    compute_metric,
    unit_sphere_area,
    panel_lonlat,
    NUM_PANELS,
)


def test_metric_symmetry_at_center():
    """面板中心 (0,0) 度量应为单位度量。"""
    m = compute_metric(0.0, 0.0)
    assert abs(m.g11 - 1.0) < 1e-15
    assert abs(m.g22 - 1.0) < 1e-15
    assert abs(m.g12 - 0.0) < 1e-15
    assert abs(m.sqrt_g - 1.0) < 1e-15


def test_metric_determinant_relation():
    """√G = sqrt(det g) 应严格成立。"""
    for xi in (-0.8, -0.3, 0.0, 0.3, 0.8):
        for eta in (-0.8, 0.2, 0.9):
            m = compute_metric(xi, eta)
            det = m.g11 * m.g22 - m.g12 * m.g12
            assert abs(m.sqrt_g - math.sqrt(det)) < 1e-14


def test_unit_sphere_area_converges_to_4pi():
    """6 面板面积积分应随分辨率收敛到 4π。"""
    # ncells=96 时相对误差应 < 3e-5（与 C++ 验证一致）
    area = unit_sphere_area(96)
    rel = abs(area - 4 * math.pi) / (4 * math.pi)
    assert rel < 3e-5, f"rel error = {rel}"


def test_area_convergence_rate():
    """面积积分随分辨率单调收敛（一阶）。"""
    a8 = abs(unit_sphere_area(8) - 4 * math.pi)
    a16 = abs(unit_sphere_area(16) - 4 * math.pi)
    a32 = abs(unit_sphere_area(32) - 4 * math.pi)
    assert a16 < a8
    assert a32 < a16


def test_panel_lonlat_north_south():
    """北/南面板中心应对应地理北极/南极。"""
    lon_n, lat_n = panel_lonlat(4, 0.0, 0.0)  # North 面板中心
    lon_s, lat_s = panel_lonlat(5, 0.0, 0.0)  # South 面板中心
    assert abs(lat_n - math.pi / 2) < 1e-12
    assert abs(lat_s + math.pi / 2) < 1e-12


def test_panel_lonlat_unit_norm():
    """面板经纬度映射应保持在单位球上（lat ∈ [-π/2, π/2]）。"""
    for panel in range(NUM_PANELS):
        for xi in (-0.9, 0.0, 0.9):
            for eta in (-0.9, 0.0, 0.9):
                lon, lat = panel_lonlat(panel, xi, eta)
                assert -math.pi / 2 <= lat <= math.pi / 2
                assert -math.pi <= lon <= math.pi
