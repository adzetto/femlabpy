import matplotlib
import numpy as np

matplotlib.use("Agg")

from matplotlib import pyplot as plt

import femlabpy as fp


def _planar_triangle_mesh():
    T = np.array([[1, 2, 3, 1]], dtype=int)
    X = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    return T, X


def test_plotelem_reduces_flat_gmsh_coordinates_to_2d():
    T, X = _planar_triangle_mesh()

    ax = fp.plotelem(T, X)

    assert ax.name != "3d"
    assert len(ax.lines) == 1
    plt.close("all")


def test_plotu_accepts_flattened_displacement_vector_with_dof():
    T, X = _planar_triangle_mesh()
    u = np.array([[0.0], [0.0], [0.1], [0.0], [0.0], [0.2]], dtype=float)

    ax = fp.plotu(T, X, u, dof=2)

    assert ax.name != "3d"
    np.testing.assert_allclose(np.asarray(ax.collections[0].get_array()), [0.1])
    plt.close("all")


def test_plotu_selects_requested_displacement_component():
    T, X = _planar_triangle_mesh()
    u = np.array([[0.0], [0.0], [0.1], [0.0], [0.0], [0.2]], dtype=float)

    ax = fp.plotu(T, X, u, dof=2, component=2)

    np.testing.assert_allclose(
        np.asarray(ax.collections[0].get_array()),
        [0.2 / 3.0],
    )
    plt.close("all")


def test_plotu_preserves_scalar_nodal_field_behavior():
    T, X = _planar_triangle_mesh()
    nodal_field = np.array([1.0, 2.0, 4.0], dtype=float)

    ax = fp.plotu(T, X[:, :2], nodal_field)

    np.testing.assert_allclose(np.asarray(ax.collections[0].get_array()), [7.0 / 3.0])
    plt.close("all")
