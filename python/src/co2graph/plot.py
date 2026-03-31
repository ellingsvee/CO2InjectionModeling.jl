"""Matplotlib-based plotting for CO2 trap graphs and simulation results."""

from __future__ import annotations

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from co2graph.graph import trap_volumes_to_grid
from co2graph.io import GraphData, LayerData


def plot_co2_heights(
    data: GraphData,
    layer_idx: int,
    volumes: np.ndarray,
    *,
    leaking: np.ndarray | None = None,
    ax: plt.Axes | None = None,
    max_co2_height: float = 20.0,
    colormap: str = "Blues",
    show_contours: bool = True,
    contour_levels: int = 10,
    major_contour_every: int = 5,
    contour_opacity: float = 0.8,
    show_extents: bool = False,
    extent_color: str = "dodgerblue",
    show_colorbar: bool = True,
    colorbar_label: str = "Column height (m)",
) -> plt.Axes:
    """Plot CO2 column heights on the spatial grid for a single layer.

    Parameters
    ----------
    data : GraphData
        Loaded export data.
    layer_idx : int
        1-indexed layer number.
    volumes : np.ndarray
        (n_traps,) CO2 volume per trap.
    leaking : np.ndarray | None
        (n_traps,) bool, for height correction on leaking traps.
    ax : plt.Axes | None
        Axes to plot on. Created if None.
    max_co2_height : float
        Colorbar maximum.
    show_extents : bool
        If True, show binary CO2 presence instead of heights.

    Returns
    -------
    plt.Axes
    """
    layer = data.layers[layer_idx - 1]
    height_grid = trap_volumes_to_grid(layer, volumes, leaking=leaking)
    return _plot_height_grid(
        data,
        layer,
        height_grid,
        ax=ax,
        max_co2_height=max_co2_height,
        colormap=colormap,
        show_contours=show_contours,
        contour_levels=contour_levels,
        major_contour_every=major_contour_every,
        contour_opacity=contour_opacity,
        show_extents=show_extents,
        extent_color=extent_color,
        show_colorbar=show_colorbar,
        colorbar_label=colorbar_label,
    )


def _plot_height_grid(
    data: GraphData,
    layer: LayerData,
    height_grid: np.ndarray,
    *,
    ax: plt.Axes | None = None,
    max_co2_height: float = 20.0,
    colormap: str = "Blues",
    show_contours: bool = True,
    contour_levels: int = 10,
    major_contour_every: int = 5,
    contour_opacity: float = 0.8,
    show_extents: bool = False,
    extent_color: str = "dodgerblue",
    show_colorbar: bool = True,
    colorbar_label: str = "Column height (m)",
) -> plt.Axes:
    """Internal: plot a height grid with optional topography contours."""
    nx, ny = layer.topography.shape
    x = np.linspace(0, (nx - 1) * data.dx, nx)
    y = np.linspace(0, (ny - 1) * data.dy, ny)
    X, Y = np.meshgrid(x, y, indexing="ij")

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))

    if show_extents:
        extent_array = (height_grid > 0).astype(float)
        cmap = mcolors.LinearSegmentedColormap.from_list(
            "extent", [(1, 1, 1, 0), mcolors.to_rgba(extent_color, 0.5)]
        )
        ax.pcolormesh(X, Y, extent_array, cmap=cmap, vmin=0, vmax=1, shading="auto")
    else:
        im = ax.pcolormesh(
            X,
            Y,
            height_grid,
            cmap=colormap,
            vmin=0,
            vmax=max_co2_height,
            shading="auto",
        )
        if show_colorbar:
            plt.colorbar(im, ax=ax, label=colorbar_label, shrink=0.8)

    if show_contours:
        topo = layer.topography
        mn, mx = topo.min(), topo.max()
        all_levels = np.linspace(mn, mx, contour_levels)
        major_levels = all_levels[1::major_contour_every]
        minor_levels = np.setdiff1d(all_levels, major_levels)

        if len(minor_levels) > 0:
            ax.contour(
                X,
                Y,
                topo,
                levels=minor_levels,
                colors="black",
                linewidths=0.5,
                alpha=contour_opacity,
            )
        if len(major_levels) > 0:
            ax.contour(
                X,
                Y,
                topo,
                levels=major_levels,
                colors="black",
                linewidths=1.5,
                alpha=contour_opacity,
            )

    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(layer.name)
    ax.set_aspect("equal")
    return ax


def plot_multi_layer(
    data: GraphData,
    volumes_per_layer: list[np.ndarray],
    *,
    leaking_per_layer: list[np.ndarray] | None = None,
    max_co2_height: float = 20.0,
    colormap: str = "Blues",
    show_contours: bool = True,
    contour_levels: int = 10,
    major_contour_every: int = 5,
    show_extents: bool = False,
    colorbar_label: str = "Column height (m)",
    figsize: tuple[int, int] | None = None,
) -> Figure:
    """Plot CO2 heights for all layers side-by-side.

    Parameters
    ----------
    data : GraphData
    volumes_per_layer : list[np.ndarray]
        Per-layer trap volumes, in layer order (deepest first).
    leaking_per_layer : list[np.ndarray] | None
        Per-layer leaking flags.
    """
    n = data.n_layers
    if figsize is None:
        figsize = (6 * n, 5)

    fig, axes = plt.subplots(1, n, figsize=figsize)
    if n == 1:
        axes = [axes]

    for i in range(n):
        layer = data.layers[i]
        leaking = leaking_per_layer[i] if leaking_per_layer else None
        height_grid = trap_volumes_to_grid(layer, volumes_per_layer[i], leaking=leaking)

        _plot_height_grid(
            data,
            layer,
            height_grid,
            ax=axes[i],
            max_co2_height=max_co2_height,
            colormap=colormap,
            show_contours=show_contours,
            contour_levels=contour_levels,
            major_contour_every=major_contour_every,
            show_extents=show_extents,
            show_colorbar=(i == n - 1),
            colorbar_label=colorbar_label,
        )
        if i > 0:
            axes[i].set_ylabel("")

    fig.tight_layout()
    return fig
