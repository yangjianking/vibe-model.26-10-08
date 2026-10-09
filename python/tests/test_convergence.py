"""平流收敛阶单元测试。"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from cubedsphere.convergence import convergence_order


def test_center2_second_order():
    """中心差分应为二阶收敛。"""
    orders = convergence_order("Center2", resolutions=(32, 64, 128))
    assert orders[-1] > 1.5, f"末阶 {orders[-1]} 应 > 1.5（二阶）"


def test_upwind3_third_order():
    """三阶上游应趋近三阶收敛。"""
    orders = convergence_order("Upwind3", resolutions=(64, 128, 256))
    assert orders[-1] > 2.5, f"末阶 {orders[-1]} 应 > 2.5（三阶）"


def test_ppm_monotone_higher_than_second():
    """PPM（单调限制）在平滑场下应优于二阶。"""
    orders = convergence_order("PPM", resolutions=(64, 128, 256))
    assert orders[-1] > 2.0, f"末阶 {orders[-1]} 应 > 2.0"
