"""标准算例与诊断单元测试。"""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from cubedsphere.test_cases import WilliamsonTC2, HeldSuarez
from cubedsphere.diagnostics import (
    global_integral,
    ConservationSnapshot,
    relative_drift,
)
from cubedsphere.config import load_yaml, get_path


def test_tc2_peak_value():
    """TC2 钟型场在中心处峰值为 1.0。"""
    tc2 = WilliamsonTC2()
    assert abs(tc2.analytic(tc2.lambda_c, tc2.phi_c, 0.0) - 1.0) < 1e-12


def test_tc2_antipode_zero():
    """TC2 在对跖点处为 0。"""
    tc2 = WilliamsonTC2()
    far = tc2.analytic(tc2.lambda_c + math.pi, tc2.phi_c, 0.0)
    assert far < 1e-12


def test_tc2_periodic_12days():
    """TC2 12 天后回到原点。"""
    tc2 = WilliamsonTC2()
    h0 = tc2.analytic(0.3, 0.2, 0.0)
    h12 = tc2.analytic(0.3, 0.2, 12.0 * 86400.0)
    assert abs(h0 - h12) < 1e-6


def test_held_suarez_equator_warmer_than_pole():
    """HS 赤道低层应暖于极地。"""
    t_eq = HeldSuarez.equilibrium_temperature(0.0, 0.0)
    t_pole = HeldSuarez.equilibrium_temperature(0.0, 1.4)
    assert t_eq > t_pole


def test_held_suarez_friction_boundary_layer():
    """HS 摩擦应在边界层最大，高层为零。"""
    kv_low = HeldSuarez.rayleigh_friction(3000.0)   # σ≈0.9
    kv_high = HeldSuarez.rayleigh_friction(24000.0)  # σ≈0.2
    assert kv_low > kv_high
    assert kv_high == 0.0


def test_global_integral_unit_field():
    """单位场全球积分 ≈ 4πR² × 总高度。"""
    ncells = 16
    nlev = 8
    radius = 6371229.0
    z_top = 30000.0
    z_half = [z_top * k / nlev for k in range(nlev + 1)]
    field = [1.0] * (6 * ncells * ncells * nlev)

    vol = global_integral(field, ncells, nlev, z_half, radius)
    area = 4 * math.pi * radius * radius
    expected = area * z_top
    rel = abs(vol - expected) / expected
    assert rel < 2e-3, f"rel error = {rel}"


def test_relative_drift_zero_for_identical():
    """相同快照的漂移应为 0。"""
    s = ConservationSnapshot(total_mass=1e18, total_theta=1e18,
                             total_energy=1e24, total_angular_momentum=1e26)
    drift = relative_drift(s, s)
    assert drift["mass"] == 0.0
    assert drift["energy"] == 0.0


def test_config_parse_simple():
    """内置 YAML 子集解析器应正确解析嵌套映射。"""
    yaml_text = """
grid:
  ncells: 48
  nhalo: 3
  radius: 6371229.0
vertical:
  nlev: 60
timeint:
  semi_implicit: true
"""
    cfg = load_yaml(yaml_text)
    assert get_path(cfg, "grid.ncells") == 48
    assert get_path(cfg, "grid.radius") == 6371229.0
    assert get_path(cfg, "vertical.nlev") == 60
    assert get_path(cfg, "timeint.semi_implicit") is True
