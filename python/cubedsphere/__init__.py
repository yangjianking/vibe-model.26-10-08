"""cubedsphere — 立方球全球大气模式 Python 后处理包。

提供立方球网格几何工具、守恒诊断、标准算例解析解与数据 I/O 辅助。

设计原则：
- 核心几何/度量计算纯 Python（标准库 math），无重依赖，可离线验证；
- 数据 I/O（NetCDF/Zarr）与绘图依赖 xarray/dask/matplotlib，lazy import，
  未安装时给出清晰报错而非导入失败。
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
