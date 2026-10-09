"""NetCDF 输出读取（xarray + Dask 惰性加载）。

依赖 xarray/dask/netCDF4，lazy import，未安装时给出清晰报错。
"""

from __future__ import annotations


class MissingDependencyError(ImportError):
    """缺少可选依赖时抛出。"""


def open_dataset(path: str, chunks: dict | None = None):
    """打开 NetCDF 输出为 xarray.Dataset（Dask 惰性加载）。

    参数：
        path    : NetCDF 文件路径
        chunks  : Dask 分块大小，如 {"time": 1, "n": "auto"}
    """
    try:
        import xarray as xr  # type: ignore
    except ImportError as e:
        raise MissingDependencyError(
            "读取 NetCDF 需要 xarray，请 `pip install xarray netCDF4 dask`"
        ) from e
    if chunks is None:
        chunks = {}
    return xr.open_dataset(path, chunks=chunks, engine="netcdf4")


def load_field(ds, name: str):
    """从数据集取出字段的 Dask 数组（惰性）。"""
    if name not in ds:
        raise KeyError(f"字段 {name!r} 不在数据集中，可用字段：{list(ds.data_vars)}")
    return ds[name]


def zonal_mean(ds, name: str) -> "xarray.DataArray":
    """计算纬向平均（惰性）。"""
    import xarray as xr  # type: ignore
    field = load_field(ds, name)
    if "lat" in field.dims:
        return field.mean(dim="panel").mean(dim="i").mean(dim="j")
    return field
