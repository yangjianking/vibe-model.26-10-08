# cubedsphere — 立方球 NWP 模式 Python 后处理包

全球大气模式的 Python 后处理工具，提供立方球网格几何、守恒诊断、标准算例
解析解与数据 I/O 辅助。

## 设计

- **核心零依赖**：网格几何（`grid.py`）、守恒诊断（`diagnostics.py`）、标准算例
  （`test_cases.py`）、配置解析（`config.py`）均为纯 Python（标准库），可离线验证。
- **可选依赖**：数据 I/O（`io.py`）依赖 xarray/dask/netCDF4；绘图（`plotting.py`）
  依赖 matplotlib。未安装时给出清晰报错，不影响核心功能。

## 安装

```bash
# 仅核心（零依赖）
pip install -e ./python

# 全功能
pip install -e "./python[all]"
```

## 快速开始

```python
from cubedsphere.grid import unit_sphere_area, compute_metric
from cubedsphere.test_cases import WilliamsonTC2

# 验证单位球面积 ≈ 4π
print(unit_sphere_area(48))  # → 12.566...

# Williamson TC2 解析解
tc2 = WilliamsonTC2()
print(tc2.analytic(lon=0.0, lat=0.0, t=0.0))
```

## 模块

| 模块 | 功能 |
|---|---|
| `grid` | gnomonic 度量、面板经纬度、面积积分 |
| `diagnostics` | 全球积分、守恒量快照、相对漂移 |
| `test_cases` | Williamson TC2 / Held–Suarez 解析解 |
| `config` | YAML 配置解析（PyYAML 优先，内置回退） |
| `io` | NetCDF 读取（xarray+Dask） |
| `plotting` | 6 面板展开图、时间序列 |

## 测试

```bash
cd python && python -m pytest -v
```
