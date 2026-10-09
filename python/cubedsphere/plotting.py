"""绘图工具（matplotlib）。

lazy import matplotlib，未安装时给出清晰报错。提供：
- 6 面板展开的全球场图（spherical projection）
- 守恒量时间序列
- 纬向平均剖面

注意：涉及地图的绘图须遵守中国地图合规要求（若绘国界/省界需用规范数据源）。
本模块仅绘制大气物理场，不涉及行政区划边界。
"""

from __future__ import annotations

from .grid import panel_lonlat, NUM_PANELS


def _plt():
    try:
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except ImportError as e:
        raise ImportError(
            "绘图需要 matplotlib，请 `pip install matplotlib`"
        ) from e


def plot_panel_field(field_flat, ncells, panel, ax=None, cmap="RdBu_r"):
    """绘制单个面板的 2D 场。

    参数：
        field_flat : 面板内扁平数组（长度 ncells×ncells，j 主序）
        ncells     : 每面板单方向单元数
        panel      : 面板索引 0–5
    """
    plt = _plt()
    import numpy as np  # type: ignore

    if ax is None:
        _, ax = plt.subplots()
    f = np.asarray(field_flat, dtype=float).reshape(ncells, ncells)
    im = ax.imshow(f, origin="lower", cmap=cmap)
    ax.set_title(f"Panel {panel}")
    ax.set_xlabel("i")
    ax.set_ylabel("j")
    return im


def plot_global_panels(field_per_panel, ncells, figsize=(12, 8), cmap="RdBu_r"):
    """6 面板展开图（2×3 子图）。"""
    plt = _plt()
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    for panel in range(NUM_PANELS):
        ax = axes[panel // 3, panel % 3]
        plot_panel_field(field_per_panel[panel], ncells, panel, ax=ax, cmap=cmap)
    fig.tight_layout()
    return fig, axes


def plot_time_series(times, values, ylabel="", title="", ax=None):
    """绘制时间序列（如守恒量漂移）。"""
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(times, values, marker="o")
    ax.set_xlabel("time [s]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return ax
