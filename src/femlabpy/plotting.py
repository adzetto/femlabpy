"""
Lightweight Matplotlib plotting helpers for quick FEM inspection.

Workflow role
-------------
These functions are not a full visualization subsystem. They are compact,
teaching-oriented helpers for checking meshes, supports, nodal loads,
deformed shapes, and contour-like stress fields directly from the FemLab-style
tables used throughout the package.

Public entry points
-------------------
- ``plotelem`` draws the undeformed mesh and optional labels.
- ``plotforces`` and ``plotbc`` overlay external loads and constraints.
- ``plotq4`` and ``plott3`` reconstruct scalar fields for quadrilateral and
  triangular meshes.
- ``plotu`` draws scalar nodal fields or displacement-derived contours.
"""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from ._helpers import as_float_array


def _axis(ax=None, ndim: int = 2):
    if ax is not None:
        return ax
    if ndim == 3:
        fig = plt.figure()
        return fig.add_subplot(111, projection="3d")
    _, axis = plt.subplots()
    return axis


def _plot_coords(X):
    """Return plotting coordinates, collapsing flat ``(x, y, z)`` meshes to 2D."""
    coords = as_float_array(X)
    if coords.ndim != 2:
        raise ValueError("X must be a two-dimensional coordinate array.")
    if coords.shape[1] == 1:
        return np.column_stack([coords[:, 0], np.zeros(coords.shape[0])])
    if coords.shape[1] > 2:
        extra = coords[:, 2:]
        if extra.size > 0 and np.allclose(extra, extra[0]):
            return coords[:, :2].copy()
    return coords


def _plot_values(u, nn: int, dof: int | None = None, component: int = 0):
    """Resolve scalar nodal values from scalar, vector, or flattened input."""
    data = as_float_array(u)
    if data.ndim == 0:
        raise ValueError("u must contain at least one nodal value.")

    flat = data.reshape(-1)
    nodal = None
    resolved_dof = dof

    if data.ndim == 1 or (data.ndim == 2 and 1 in data.shape):
        if dof is None:
            if flat.size != nn:
                raise ValueError(
                    "When dof is omitted, u must contain one scalar value per node."
                )
            return flat
        if flat.size == nn:
            return flat
        if flat.size != nn * dof:
            raise ValueError(
                f"u must contain either {nn} nodal values or {nn * dof} "
                f"flattened values for dof={dof}."
            )
        nodal = flat.reshape(nn, dof)
        resolved_dof = dof
    elif data.ndim == 2:
        if data.shape[0] == nn:
            nodal = data
        elif data.shape[1] == nn:
            nodal = data.T
        else:
            raise ValueError(
                "u must be a nodal array with shape (nn,), (nn, ncomp), or a "
                "flattened vector with length nn * dof."
            )
        resolved_dof = nodal.shape[1] if dof is None else dof
        if nodal.shape[1] < resolved_dof:
            raise ValueError(
                f"u provides {nodal.shape[1]} nodal components, but dof={resolved_dof} "
                "was requested."
            )
    else:
        raise ValueError("u must be one- or two-dimensional.")

    if nodal is None or resolved_dof is None:
        raise ValueError("Could not resolve nodal values for plotting.")
    if component < 0:
        raise ValueError("component must be non-negative.")
    if nodal.shape[1] == 1:
        return nodal[:, 0]
    if component == 0:
        return np.linalg.norm(nodal[:, :resolved_dof], axis=1)
    if component > resolved_dof:
        raise ValueError(
            f"component={component} is out of range for dof={resolved_dof}."
        )
    return nodal[:, component - 1]


def plotelem(
    T, X, line_style: str = "k-", nonum: bool = False, noelem: bool = False, ax=None
):
    """Plot the undeformed mesh and optionally annotate node or element numbers.

    Algorithm
    ---------
    1. Normalize the plotting coordinates, reducing flat Gmsh ``(x, y, z)``
       input to a planar ``(x, y)`` view.
    2. Loop over the element rows and draw the corresponding polygon or line.
    3. Optionally annotate node ids or element ids.
    """
    topology = as_float_array(T).astype(int)
    coords = _plot_coords(X)
    ndim = coords.shape[1]
    ax = _axis(ax, ndim)
    for idx, row in enumerate(topology, start=1):
        nodes = row[:-1] - 1
        order = list(nodes)
        if len(order) >= 3:
            order.append(order[0])
        points = coords[order]
        if ndim == 2:
            ax.plot(points[:, 0], points[:, 1], line_style)
            if nonum:
                for node in nodes:
                    ax.text(coords[node, 0], coords[node, 1], str(node + 1), fontsize=8)
            if noelem:
                center = coords[nodes].mean(axis=0)
                ax.text(center[0], center[1], str(idx), fontsize=8)
        else:
            ax.plot(points[:, 0], points[:, 1], points[:, 2], line_style)
    if ndim == 2:
        ax.set_aspect("equal", adjustable="box")
    return ax


def plotforces(T, X, P, ax=None):
    """Plot nodal loads as arrows on a 2D mesh view.

    Algorithm
    ---------
    1. Extract element nodal coordinates.
    2. Build a `matplotlib.collections.PolyCollection`.
    3. Map internal stresses to a colormap using `vmin` and `vmax`.
    """
    _ = T
    coords = _plot_coords(X)
    loads = as_float_array(P)
    ax = _axis(ax, 2)
    if loads.size == 0:
        return ax
    max_force = max(np.max(np.abs(loads[:, 1:3])), 1.0)
    span = np.ptp(coords[:, :2], axis=0)
    scale = 0.5 * max(span.max(), 1.0)
    for row in loads:
        node = int(row[0]) - 1
        x0, y0 = coords[node, :2]
        dx = scale * row[1] / max_force
        dy = scale * row[2] / max_force
        ax.arrow(
            x0, y0, dx, dy, color="tab:green", width=0.002, length_includes_head=True
        )
    return ax


def plotbc(T, X, C, ax=None):
    """Plot prescribed boundary conditions on a 2D mesh view.

    Algorithm
    ---------
    1. Extract element nodal coordinates.
    2. Build a `matplotlib.collections.PolyCollection`.
    3. Map internal stresses to a colormap using `vmin` and `vmax`.
    """
    _ = T
    coords = _plot_coords(X)
    constraints = as_float_array(C)
    ax = _axis(ax, 2)
    span = np.ptp(coords[:, :2], axis=0)
    scale = 0.1 * max(span.max(), 1.0)
    for row in constraints:
        node = int(row[0]) - 1
        dof = int(row[1])
        x0, y0 = coords[node, :2]
        value = row[-1]
        if value == 0:
            if dof == 1:
                ax.scatter([x0 - scale], [y0], color="black", marker="s", s=20)
            else:
                ax.scatter([x0], [y0 - scale], color="black", marker="s", s=20)
        else:
            dx = value if dof == 1 else 0.0
            dy = value if dof == 2 else 0.0
            ax.arrow(
                x0, y0, dx, dy, color="tab:red", width=0.002, length_includes_head=True
            )
    return ax


def _triangulate_quads(T):
    topology = as_float_array(T).astype(int)
    triangles = []
    for row in topology:
        nodes = row[:-1]
        triangles.append([nodes[0], nodes[1], nodes[2]])
        triangles.append([nodes[0], nodes[2], nodes[3]])
    return np.asarray(triangles, dtype=int)


def plotq4(T, X, S, scomp: int, ax=None):
    """Plot a contour field reconstructed from Q4 Gauss-point results.

    Algorithm
    ---------
    1. Extract element nodal coordinates.
    2. Build a `matplotlib.collections.PolyCollection`.
    3. Map internal stresses to a colormap using `vmin` and `vmax`.
    """
    topology = as_float_array(T).astype(int)
    coords = _plot_coords(X)
    values = as_float_array(S)
    ncomp = values.shape[1] // 4
    if scomp > ncomp:
        raise ValueError(f"Requested component {scomp} is not available.")
    r = np.array([-1.0, 1.0], dtype=float) * np.sqrt(3.0)
    N = np.zeros((4, 4), dtype=float)
    for i in range(2):
        for j in range(2):
            gp = i + 3 * j - 2 * i * j
            N[gp] = (
                np.array(
                    [
                        (1.0 - r[i]) * (1.0 - r[j]),
                        (1.0 + r[i]) * (1.0 - r[j]),
                        (1.0 + r[i]) * (1.0 + r[j]),
                        (1.0 - r[i]) * (1.0 + r[j]),
                    ],
                    dtype=float,
                )
                / 4.0
            )
    nodal_values = np.zeros((coords.shape[0],), dtype=float)
    counts = np.zeros((coords.shape[0],), dtype=float)
    component_index = scomp - 1
    for element, row in enumerate(topology):
        nodes = row[:-1] - 1
        gauss_values = values[element, component_index::ncomp]
        node_values = N @ gauss_values.reshape(-1, 1)
        nodal_values[nodes] += node_values.ravel()
        counts[nodes] += 1.0
    nodal_values /= np.maximum(counts, 1.0)
    triangles = _triangulate_quads(topology) - 1
    ax = _axis(ax, 2)
    trip = ax.tripcolor(
        coords[:, 0], coords[:, 1], triangles, nodal_values, shading="gouraud"
    )
    ax.set_aspect("equal", adjustable="box")
    plt.colorbar(trip, ax=ax)
    return ax


def plott3(T, X, S, scomp: int, ax=None):
    """Plot a contour field from T3 element results.

    Algorithm
    ---------
    1. Extract element nodal coordinates.
    2. Build a `matplotlib.collections.PolyCollection`.
    3. Map internal stresses to a colormap using `vmin` and `vmax`.
    """
    topology = as_float_array(T).astype(int)
    coords = _plot_coords(X)
    values = as_float_array(S)
    if scomp > values.shape[1]:
        raise ValueError(f"Requested component {scomp} is not available.")
    nodal_values = np.zeros((coords.shape[0],), dtype=float)
    counts = np.zeros((coords.shape[0],), dtype=float)
    for element, row in enumerate(topology):
        nodes = row[:-1] - 1
        nodal_values[nodes] += values[element, scomp - 1]
        counts[nodes] += 1.0
    nodal_values /= np.maximum(counts, 1.0)
    triangles = topology[:, :-1] - 1
    ax = _axis(ax, 2)
    trip = ax.tripcolor(
        coords[:, 0], coords[:, 1], triangles, nodal_values, shading="flat"
    )
    ax.set_aspect("equal", adjustable="box")
    plt.colorbar(trip, ax=ax)
    return ax


def plotu(T, X, u, dof: int | None = None, component: int = 0, ax=None):
    """Plot a scalar nodal field or displacement-derived contour.

    Parameters
    ----------
    T : array_like
        Element topology table with one-based node ids and a trailing property id.
    X : array_like
        Nodal coordinates. Flat Gmsh ``(x, y, z)`` meshes are plotted in 2D.
    u : array_like
        Either one scalar value per node, a ``(nn, ncomp)`` nodal field, or a
        flattened global vector with length ``nn * dof``.
    dof : int, optional
        Degrees of freedom per node when ``u`` is passed as a flattened global
        vector. When omitted, ``u`` is treated as a scalar nodal field unless it
        is already shaped as ``(nn, ncomp)``.
    component : int, default 0
        Which nodal component to plot when ``u`` contains multiple values per
        node. Use ``0`` for magnitude, ``1`` for the first component,
        ``2`` for the second component, and so on.
    ax : matplotlib.axes.Axes, optional
        Existing Matplotlib axes. A new axes object is created when omitted.

    Algorithm
    ---------
    1. Normalize the coordinates for plotting, collapsing flat ``(x, y, z)``
       geometry to a planar view when appropriate.
    2. Resolve one scalar value per node from ``u``. Flattened vectors use
       ``dof`` and ``component`` to extract a component or displacement magnitude.
    3. Color each element by the mean of its nodal values.
    """
    topology = as_float_array(T).astype(int)
    coords = _plot_coords(X)
    values = _plot_values(u, coords.shape[0], dof=dof, component=component)
    ndim = coords.shape[1]
    ax = _axis(ax, ndim)
    if ndim == 2:
        polygons = []
        colors = []
        for row in topology:
            nodes = row[:-1] - 1
            polygons.append(coords[nodes, :2])
            colors.append(values[nodes].mean())
        collection = PolyCollection(
            polygons, array=np.asarray(colors), cmap="viridis", edgecolors="k"
        )
        ax.add_collection(collection)
        ax.autoscale()
        ax.set_aspect("equal", adjustable="box")
        plt.colorbar(collection, ax=ax)
        return ax

    polygons3d = []
    colors = []
    for row in topology:
        nodes = row[:-1] - 1
        polygons3d.append(coords[nodes, :3])
        colors.append(values[nodes].mean())
    collection = Poly3DCollection(
        polygons3d, array=np.asarray(colors), cmap="viridis", edgecolors="k"
    )
    ax.add_collection3d(collection)
    ax.auto_scale_xyz(coords[:, 0], coords[:, 1], coords[:, 2])
    plt.colorbar(collection, ax=ax)
    return ax


__all__ = ["plotbc", "plotelem", "plotforces", "plotq4", "plott3", "plotu"]
