"""
Three-dimensional solid element kernels for tetrahedra and hexahedra.

Workflow role
-------------
This module contains the 3D structural kernels used when the model carries
three translational degrees of freedom per node. It covers both simplex and
brick topologies and includes stiffness, internal-force, and mass-matrix
construction at element and assembled levels.

Public entry points
-------------------
- ``keT4e`` and ``qeT4e`` operate on one tetrahedral element.
- ``keh8e`` and ``qeh8e`` operate on one hexahedral brick.
- ``kT4e``, ``qT4e``, ``kh8e``, and ``qh8e`` assemble those responses globally.
- ``meT4e``, ``mT4e``, ``meh8e``, and ``mh8e`` provide the dynamic mass terms.
"""

from __future__ import annotations

import numpy as np

from .._helpers import as_float_array, element_dof_indices, is_sparse

try:
    import scipy.sparse as sp
except ImportError:  # pragma: no cover
    sp = None


def _elastic3d_matrix(Ge):
    """
    Build the isotropic 3D constitutive matrix for one material row.

    Parameters
    ----------
    Ge:
        Material row containing at least ``[E, nu]``.

    Returns
    -------
    ndarray
        ``6 x 6`` 3D elastic matrix in Voigt order.
    """
    props = as_float_array(Ge).reshape(-1)
    E = props[0]
    nu = props[1]
    return (
        E
        / ((1.0 + nu) * (1.0 - 2.0 * nu))
        * np.array(
            [
                [1.0 - nu, nu, nu, 0.0, 0.0, 0.0],
                [nu, 1.0 - nu, nu, 0.0, 0.0, 0.0],
                [nu, nu, 1.0 - nu, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, (1.0 - 2.0 * nu) / 2.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, (1.0 - 2.0 * nu) / 2.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, (1.0 - 2.0 * nu) / 2.0],
            ],
            dtype=float,
        )
    )


def _elastic3d_matrix_batch(Ge):
    """
    Vectorize :func:`_elastic3d_matrix` over multiple materials.

    Parameters
    ----------
    Ge:
        Material table with one row per element.

    Returns
    -------
    ndarray
        Batched constitutive matrices with shape ``(n, 6, 6)``.
    """
    props = as_float_array(Ge)
    if props.ndim == 1:
        props = props.reshape(1, -1)
    E = props[:, 0]
    nu = props[:, 1]
    D = np.zeros((props.shape[0], 6, 6), dtype=float)
    D[:, 0, 0] = 1.0 - nu
    D[:, 0, 1] = nu
    D[:, 0, 2] = nu
    D[:, 1, 0] = nu
    D[:, 1, 1] = 1.0 - nu
    D[:, 1, 2] = nu
    D[:, 2, 0] = nu
    D[:, 2, 1] = nu
    D[:, 2, 2] = 1.0 - nu
    D[:, 3, 3] = (1.0 - 2.0 * nu) / 2.0
    D[:, 4, 4] = (1.0 - 2.0 * nu) / 2.0
    D[:, 5, 5] = (1.0 - 2.0 * nu) / 2.0
    return D * (E / ((1.0 + nu) * (1.0 - 2.0 * nu)))[:, None, None]


def _solid_B(dN):
    """
    Assemble the 3D strain-displacement matrix for one solid element.

    Parameters
    ----------
    dN:
        Global shape-function gradients with shape ``(3, nnodes)``.

    Returns
    -------
    ndarray
        ``6 x (3 * nnodes)`` strain-displacement matrix.
    """
    nnodes = dN.shape[1]
    B = np.zeros((6, nnodes * 3), dtype=float)
    B[0, 0::3] = dN[0]
    B[1, 1::3] = dN[1]
    B[2, 2::3] = dN[2]
    B[3, 0::3] = dN[1]
    B[3, 1::3] = dN[0]
    B[4, 1::3] = dN[2]
    B[4, 2::3] = dN[1]
    B[5, 0::3] = dN[2]
    B[5, 2::3] = dN[0]
    return B


def _solid_B_batch(dN):
    """
    Vectorize :func:`_solid_B` over one or more integration points.

    Parameters
    ----------
    dN:
        Batched global shape-function gradients.

    Returns
    -------
    ndarray
        Batched strain-displacement matrices.
    """
    nnodes = dN.shape[-1]
    B = np.zeros(dN.shape[:-2] + (6, nnodes * 3), dtype=float)
    B[..., 0, 0::3] = dN[..., 0, :]
    B[..., 1, 1::3] = dN[..., 1, :]
    B[..., 2, 2::3] = dN[..., 2, :]
    B[..., 3, 0::3] = dN[..., 1, :]
    B[..., 3, 1::3] = dN[..., 0, :]
    B[..., 4, 1::3] = dN[..., 2, :]
    B[..., 4, 2::3] = dN[..., 1, :]
    B[..., 5, 0::3] = dN[..., 2, :]
    B[..., 5, 2::3] = dN[..., 0, :]
    return B


def keT4e(Xe, Ge):
    """
    Compute the stiffness matrix for a 4-node tetrahedral solid element.

    Parameters
    ----------
    Xe:
        Tetrahedral nodal coordinates with shape ``(4, 3)``.
    Ge:
        Material row ``[E, nu]``.

    Returns
    -------
    ndarray
        ``12 x 12`` element stiffness matrix.

    Mathematical Formulation
    ------------------------
    For a linear tetrahedron, the strain-displacement matrix `B` is constant
    within the element, so the stiffness matrix is obtained from the exact
    volume integral of `B.T @ C @ B`.

    Algorithm
    ---------
    1. Extract nodal coordinates and calculate the Jacobian matrix.
    2. Compute the shape-function derivatives and the `6 x 12` matrix `B`.
    3. Construct the `6 x 6` constitutive matrix `C`.
    4. Compute the element stiffness matrix.
    """
    Xe = as_float_array(Xe)
    dN = np.array(
        [[1.0, 0.0, 0.0, -1.0], [0.0, 1.0, 0.0, -1.0], [0.0, 0.0, 1.0, -1.0]],
        dtype=float,
    )
    J = dN @ Xe
    dN = np.linalg.solve(J, dN)
    B = _solid_B(dN)
    D = _elastic3d_matrix(Ge)
    return 2.0 * (B.T @ D @ B) * np.linalg.det(J)


def qeT4e(Xe, Ge, Ue):
    """
    Recover stress and strain for one tetrahedral solid element.

    Parameters
    ----------
    Xe:
        Tetrahedral nodal coordinates with shape ``(4, 3)``.
    Ge:
        Material row ``[E, nu]``.
    Ue:
        Element displacement vector.

    Returns
    -------
    tuple[ndarray, ndarray, ndarray]
        Element internal-force vector, stress vector, and strain vector.

    Mathematical Formulation
    ------------------------
    Strain is computed as `Ee = B @ Ue`, stress is `Se = C @ Ee`, and the
    internal-force vector is obtained from the element volume times
    `B.T @ Se`.

    Algorithm
    ---------
    1. Compute the `B` matrix and the constitutive matrix `C`.
    2. Evaluate the element strain and stress vectors.
    3. Compute the element internal-force vector.
    """
    Xe = as_float_array(Xe)
    Ue = as_float_array(Ue).reshape(-1, 1)
    dN = np.array(
        [[1.0, 0.0, 0.0, -1.0], [0.0, 1.0, 0.0, -1.0], [0.0, 0.0, 1.0, -1.0]],
        dtype=float,
    )
    J = dN @ Xe
    dN = np.linalg.solve(J, dN)
    B = _solid_B(dN)
    D = _elastic3d_matrix(Ge)
    Ee = (B @ Ue).reshape(-1)
    Se = Ee @ D
    qe = (B.T @ Se.reshape(-1, 1)) * np.linalg.det(J)
    return qe, Se, Ee


def kT4e(K, T, X, G):
    """
    Assemble T4 solid element stiffness contributions into the global matrix.

    Parameters
    ----------
    K:
        Global stiffness matrix.
    T:
        Topology table ``[n1, n2, n3, n4, mat_id]``.
    X:
        Nodal coordinates.
    G:
        Material table.

    Returns
    -------
    ndarray or sparse matrix
        Updated global stiffness matrix.

    Mathematical Formulation
    ------------------------
    The global stiffness matrix is assembled from the element matrices in the
    standard form `L_e.T @ Ke @ L_e`, where `L_e` is the Boolean connectivity
    operator for element `e`.

    Algorithm
    ---------
    1. Loop over all elements or compute in a vectorized batch.
    2. Compute element stiffness matrices using batched operations.
    3. Determine global degrees of freedom indices for each element.
    4. Scatter the element stiffness entries into the global stiffness matrix.
    """
    topology = as_float_array(T)
    coordinates = as_float_array(X)
    nodes = topology[:, :4].astype(int) - 1
    materials = as_float_array(G)[topology[:, -1].astype(int) - 1]
    Xe = coordinates[nodes]
    dN = np.array(
        [[1.0, 0.0, 0.0, -1.0], [0.0, 1.0, 0.0, -1.0], [0.0, 0.0, 1.0, -1.0]],
        dtype=float,
    )
    J = np.einsum("ij,ejk->eik", dN, Xe)
    dN_global = np.linalg.solve(J, np.broadcast_to(dN, (topology.shape[0],) + dN.shape))
    B = _solid_B_batch(dN_global)
    D = _elastic3d_matrix_batch(materials)
    detJ = np.linalg.det(J)
    element_matrices = (
        2.0
        * detJ[:, None, None]
        * np.einsum("eik,ekl,elj->eij", B.transpose(0, 2, 1), D, B)
    )
    indices = element_dof_indices(nodes, 3, one_based=False)
    if is_sparse(K) and sp is not None:
        scatter_rows = np.broadcast_to(
            indices[:, :, None], element_matrices.shape
        ).reshape(-1)
        scatter_cols = np.broadcast_to(
            indices[:, None, :], element_matrices.shape
        ).reshape(-1)
        delta = sp.coo_matrix(
            (element_matrices.reshape(-1), (scatter_rows, scatter_cols)),
            shape=K.shape,
            dtype=float,
        )
        return (K.tocsr() + delta.tocsr()).tolil()
    np.add.at(K, (indices[:, :, None], indices[:, None, :]), element_matrices)
    return K


def qT4e(q, T, X, G, u):
    """
    Recover T4 solid stresses and assemble internal forces.

    Parameters
    ----------
    q:
        Global internal-force vector.
    T:
        Topology table ``[n1, n2, n3, n4, mat_id]``.
    X:
        Nodal coordinates.
    G:
        Material table.
    u:
        Global displacement vector.

    Returns
    -------
    tuple[ndarray, ndarray, ndarray]
        Updated internal-force vector, element stresses, and element strains.

    Mathematical Formulation
    ------------------------
    The global internal-force vector is assembled from the element vectors in
    the standard form `L_e.T @ qe`, with strains and stresses evaluated
    element by element.

    Algorithm
    ---------
    1. Extract element displacements from the global displacement vector.
    2. Compute batched `B` matrices and element strains.
    3. Compute batched stresses from the elastic constitutive matrices.
    4. Scatter element internal forces into the global force vector.
    """
    topology = as_float_array(T)
    coordinates = as_float_array(X)
    U = as_float_array(u).reshape(coordinates.shape[0], coordinates.shape[1])
    nodes = topology[:, :4].astype(int) - 1
    materials = as_float_array(G)[topology[:, -1].astype(int) - 1]
    Xe = coordinates[nodes]
    dN = np.array(
        [[1.0, 0.0, 0.0, -1.0], [0.0, 1.0, 0.0, -1.0], [0.0, 0.0, 1.0, -1.0]],
        dtype=float,
    )
    J = np.einsum("ij,ejk->eik", dN, Xe)
    dN_global = np.linalg.solve(J, np.broadcast_to(dN, (topology.shape[0],) + dN.shape))
    B = _solid_B_batch(dN_global)
    D = _elastic3d_matrix_batch(materials)
    detJ = np.linalg.det(J)
    element_displacements = U[nodes].reshape(topology.shape[0], -1)
    E = np.einsum("eij,ej->ei", B, element_displacements)
    S = np.einsum("ei,eij->ej", E, D)
    element_vectors = detJ[:, None] * np.einsum("eij,ej->ei", B.transpose(0, 2, 1), S)
    indices = element_dof_indices(nodes, coordinates.shape[1], one_based=False)
    np.add.at(q[:, 0], indices.reshape(-1), element_vectors.reshape(-1))
    return q, S, E


def _hexa_dN(r_i: float, r_j: float, r_k: float):
    """
    Evaluate H8 parent-space shape-function derivatives.

    Parameters
    ----------
    r_i, r_j, r_k:
        Parent-space Gauss-point coordinates ``(xi, eta, zeta)``.

    Returns
    -------
    ndarray
        Parent-space derivatives with shape ``(3, 8)``.
    """
    dNi = (
        np.array(
            [
                -(1.0 - r_j) * (1.0 - r_k),
                (1.0 - r_j) * (1.0 - r_k),
                (1.0 + r_j) * (1.0 - r_k),
                -(1.0 + r_j) * (1.0 - r_k),
                -(1.0 - r_j) * (1.0 + r_k),
                (1.0 - r_j) * (1.0 + r_k),
                (1.0 + r_j) * (1.0 + r_k),
                -(1.0 + r_j) * (1.0 + r_k),
            ],
            dtype=float,
        )
        / 8.0
    )
    dNj = (
        np.array(
            [
                -(1.0 - r_i) * (1.0 - r_k),
                -(1.0 + r_i) * (1.0 - r_k),
                (1.0 + r_i) * (1.0 - r_k),
                (1.0 - r_i) * (1.0 - r_k),
                -(1.0 - r_i) * (1.0 + r_k),
                -(1.0 + r_i) * (1.0 + r_k),
                (1.0 + r_i) * (1.0 + r_k),
                (1.0 - r_i) * (1.0 + r_k),
            ],
            dtype=float,
        )
        / 8.0
    )
    dNk = (
        np.array(
            [
                -(1.0 - r_i) * (1.0 - r_j),
                -(1.0 + r_i) * (1.0 - r_j),
                -(1.0 + r_i) * (1.0 + r_j),
                -(1.0 - r_i) * (1.0 + r_j),
                (1.0 - r_i) * (1.0 - r_j),
                (1.0 + r_i) * (1.0 - r_j),
                (1.0 + r_i) * (1.0 + r_j),
                (1.0 - r_i) * (1.0 + r_j),
            ],
            dtype=float,
        )
        / 8.0
    )
    return np.vstack([dNi, dNj, dNk])


def _hexa_dN_batch(points):
    """
    Vectorize :func:`_hexa_dN` over a batch of parent-space points.

    Parameters
    ----------
    points:
        Parent-space coordinates with shape ``(n, 3)``.

    Returns
    -------
    ndarray
        Batched derivatives with shape ``(n, 3, 8)``.
    """
    points = as_float_array(points)
    r_i = points[:, 0][:, None]
    r_j = points[:, 1][:, None]
    r_k = points[:, 2][:, None]
    dNi = (
        np.hstack(
            [
                -(1.0 - r_j) * (1.0 - r_k),
                (1.0 - r_j) * (1.0 - r_k),
                (1.0 + r_j) * (1.0 - r_k),
                -(1.0 + r_j) * (1.0 - r_k),
                -(1.0 - r_j) * (1.0 + r_k),
                (1.0 - r_j) * (1.0 + r_k),
                (1.0 + r_j) * (1.0 + r_k),
                -(1.0 + r_j) * (1.0 + r_k),
            ]
        )
        / 8.0
    )
    dNj = (
        np.hstack(
            [
                -(1.0 - r_i) * (1.0 - r_k),
                -(1.0 + r_i) * (1.0 - r_k),
                (1.0 + r_i) * (1.0 - r_k),
                (1.0 - r_i) * (1.0 - r_k),
                -(1.0 - r_i) * (1.0 + r_k),
                -(1.0 + r_i) * (1.0 + r_k),
                (1.0 + r_i) * (1.0 + r_k),
                (1.0 - r_i) * (1.0 + r_k),
            ]
        )
        / 8.0
    )
    dNk = (
        np.hstack(
            [
                -(1.0 - r_i) * (1.0 - r_j),
                -(1.0 + r_i) * (1.0 - r_j),
                -(1.0 + r_i) * (1.0 + r_j),
                -(1.0 - r_i) * (1.0 + r_j),
                (1.0 - r_i) * (1.0 - r_j),
                (1.0 + r_i) * (1.0 - r_j),
                (1.0 + r_i) * (1.0 + r_j),
                (1.0 - r_i) * (1.0 + r_j),
            ]
        )
        / 8.0
    )
    return np.stack([dNi, dNj, dNk], axis=1)


def keh8e(Xe, Ge):
    """
    Compute the stiffness matrix for an 8-node hexahedral solid element.

    Parameters
    ----------
    Xe:
        Hexahedral nodal coordinates with shape ``(8, 3)``.
    Ge:
        Material row ``[E, nu]``.

    Returns
    -------
    ndarray
        ``24 x 24`` element stiffness matrix integrated with ``2 x 2 x 2``
        Gauss quadrature.

    Mathematical Formulation
    ------------------------
    The stiffness matrix for an H8 element uses `2 x 2 x 2` Gauss integration.
    At each of the eight Gauss points, the routine evaluates
    `B_i.T @ C @ B_i * det(J_i)` and sums the contributions.

    Algorithm
    ---------
    1. Loop over 8 Gauss points.
    2. Compute the `3 x 3` Jacobian matrix.
    3. Compute the `6 x 24` matrix `B`.
    4. Accumulate the stiffness matrix contributions using Einstein summation.
    """
    Xe = as_float_array(Xe)
    if Xe.shape[0] == 20:
        raise NotImplementedError(
            "20-node hexahedra are not supported yet, matching the Scilab limitation."
        )
    D = _elastic3d_matrix(Ge)
    gauss_points = np.array(
        [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ],
        dtype=float,
    ) / np.sqrt(3.0)
    dN = _hexa_dN_batch(gauss_points)
    Jt = np.einsum("gik,kj->gij", dN, Xe)
    dN_global = np.linalg.solve(Jt, dN)
    B = _solid_B_batch(dN_global)
    return np.einsum(
        "g,gik,kl,glj->ij",
        np.linalg.det(Jt),
        B.transpose(0, 2, 1),
        D,
        B,
    )


def qeh8e(Xe, Ge, Ue):
    """
    Recover stress and strain at the eight H8 Gauss points.

    Parameters
    ----------
    Xe:
        Hexahedral nodal coordinates with shape ``(8, 3)``.
    Ge:
        Material row ``[E, nu]``.
    Ue:
        Element displacement vector.

    Returns
    -------
    tuple[ndarray, ndarray, ndarray]
        Element internal-force vector, Gauss-point stresses, and Gauss-point
        strains.

    Mathematical Formulation
    ------------------------
    Strains and stresses are evaluated at the `2 x 2 x 2` Gauss points using
    `epsilon = B @ ue` and `sigma = C @ epsilon`.

    Algorithm
    ---------
    1. Loop over 8 Gauss points.
    2. Compute the `3 x 3` Jacobian matrix and the `B` matrix at each point.
    3. Evaluate strains and stresses at each integration point.
    4. Integrate the internal forces over the volume.
    """
    Xe = as_float_array(Xe)
    if Xe.shape[0] == 20:
        raise NotImplementedError(
            "20-node hexahedra are not supported yet, matching the Scilab limitation."
        )
    D = _elastic3d_matrix(Ge)
    Ue = as_float_array(Ue).reshape(-1, 1)
    gauss_points = np.array(
        [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ],
        dtype=float,
    ) / np.sqrt(3.0)
    dN = _hexa_dN_batch(gauss_points)
    Jt = np.einsum("gik,kj->gij", dN, Xe)
    dN_global = np.linalg.solve(Jt, dN)
    B = _solid_B_batch(dN_global)
    Ee = np.einsum("gij,jk->gi", B, Ue).reshape(8, 6)
    Se = np.einsum("gi,ij->gj", Ee, D)
    qe = np.einsum(
        "g,gij,gj->i",
        np.linalg.det(Jt),
        B.transpose(0, 2, 1),
        Se,
    ).reshape(-1, 1)
    return qe, Se, Ee


def kh8e(K, T, X, G):
    """
    Assemble H8 solid element stiffness contributions into the global matrix.

    Parameters
    ----------
    K:
        Global stiffness matrix.
    T:
        Topology table ``[n1, ..., n8, mat_id]``.
    X:
        Nodal coordinates.
    G:
        Material table.

    Returns
    -------
    ndarray or sparse matrix
        Updated global stiffness matrix.

    Mathematical Formulation
    ------------------------
    Global stiffness is formed by summing the element contributions in the
    usual form `L_e.T @ Ke @ L_e`, using `2 x 2 x 2` numerical integration for
    the hexahedral elements.

    Algorithm
    ---------
    1. Loop over all elements or perform batched computations.
    2. Accumulate integrations over 8 Gauss points per element.
    3. Form batched `B` matrices and compute the stiffness components.
    4. Scatter the element stiffness matrices into the global matrix.
    """
    topology = as_float_array(T)
    coordinates = as_float_array(X)
    nodes = topology[:, :8].astype(int) - 1
    materials = as_float_array(G)[topology[:, -1].astype(int) - 1]
    Xe = coordinates[nodes]
    gauss_points = np.array(
        [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ],
        dtype=float,
    ) / np.sqrt(3.0)
    dN = _hexa_dN_batch(gauss_points)
    Jt = np.einsum("gik,ekj->egij", dN, Xe)
    dN_global = np.linalg.solve(
        Jt, np.broadcast_to(dN, (topology.shape[0],) + dN.shape)
    )
    B = _solid_B_batch(dN_global)
    D = _elastic3d_matrix_batch(materials)
    detJ = np.linalg.det(Jt)
    element_matrices = np.einsum(
        "eg,egik,ekl,eglj->eij",
        detJ,
        B.transpose(0, 1, 3, 2),
        D,
        B,
    )
    indices = element_dof_indices(nodes, 3, one_based=False)
    if is_sparse(K) and sp is not None:
        scatter_rows = np.broadcast_to(
            indices[:, :, None], element_matrices.shape
        ).reshape(-1)
        scatter_cols = np.broadcast_to(
            indices[:, None, :], element_matrices.shape
        ).reshape(-1)
        delta = sp.coo_matrix(
            (element_matrices.reshape(-1), (scatter_rows, scatter_cols)),
            shape=K.shape,
            dtype=float,
        )
        return (K.tocsr() + delta.tocsr()).tolil()
    np.add.at(K, (indices[:, :, None], indices[:, None, :]), element_matrices)
    return K


def qh8e(q, T, X, G, u):
    """
    Recover H8 solid stresses and assemble internal forces.

    Parameters
    ----------
    q:
        Global internal-force vector.
    T:
        Topology table ``[n1, ..., n8, mat_id]``.
    X:
        Nodal coordinates.
    G:
        Material table.
    u:
        Global displacement vector.

    Returns
    -------
    tuple[ndarray, ndarray, ndarray]
        Updated internal-force vector, flattened Gauss-point stresses, and
        flattened Gauss-point strains.

    Mathematical Formulation
    ------------------------
    The global internal-force vector is assembled from the element vectors,
    where each element vector is the sum over the eight Gauss points of
    `B_i.T @ sigma_i * det(J_i)`.

    Algorithm
    ---------
    1. Retrieve global displacements and evaluate localized counterparts.
    2. Form the `B` matrices using the `2 x 2 x 2` Gauss points.
    3. Map displacements to element strains and element stresses.
    4. Assemble localized internal element forces into the global vector.
    """
    topology = as_float_array(T)
    coordinates = as_float_array(X)
    U = as_float_array(u).reshape(coordinates.shape[0], coordinates.shape[1])
    nodes = topology[:, :8].astype(int) - 1
    materials = as_float_array(G)[topology[:, -1].astype(int) - 1]
    Xe = coordinates[nodes]
    gauss_points = np.array(
        [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ],
        dtype=float,
    ) / np.sqrt(3.0)
    dN = _hexa_dN_batch(gauss_points)
    Jt = np.einsum("gik,ekj->egij", dN, Xe)
    dN_global = np.linalg.solve(
        Jt, np.broadcast_to(dN, (topology.shape[0],) + dN.shape)
    )
    B = _solid_B_batch(dN_global)
    D = _elastic3d_matrix_batch(materials)
    detJ = np.linalg.det(Jt)
    element_displacements = U[nodes].reshape(topology.shape[0], -1)
    E = np.einsum("egij,ej->egi", B, element_displacements).reshape(
        topology.shape[0], -1
    )
    strain_gp = E.reshape(topology.shape[0], 8, 6)
    stress_gp = np.einsum("egi,eij->egj", strain_gp, D)
    element_vectors = np.einsum(
        "eg,egij,egj->ei",
        detJ,
        B.transpose(0, 1, 3, 2),
        stress_gp,
    )
    indices = element_dof_indices(nodes, coordinates.shape[1], one_based=False)
    np.add.at(q[:, 0], indices.reshape(-1), element_vectors.reshape(-1))
    return q, stress_gp.reshape(topology.shape[0], -1), E


def meT4e(Xe, Ge, *, lumped: bool = False):
    """
    Compute the element mass matrix for a 4-node tetrahedral solid element.

    Consistent (12x12) — analytical formula:
        M = (rho * V / 20) * [[2I, I, I, I],
                               [I, 2I, I, I],
                               [I, I, 2I, I],
                               [I, I, I, 2I]]

    Lumped (diagonal):
        M = (rho * V / 4) * I_12

    Parameters
    ----------
    Xe : array_like, shape (4, 3)
        Nodal coordinates.
    Ge : array_like
        Material row ``[E, nu, rho]``.  If ``rho`` is omitted, defaults to 1.
    lumped : bool
        If True, return diagonally lumped mass.

    Returns
    -------
    Me : ndarray, shape (12, 12)

    Mathematical Formulation
    ------------------------
    The consistent mass matrix integrates `rho * N.T @ N` exactly over the
    tetrahedral volume. The lumped version distributes the total mass equally
    among the four nodes.

    Algorithm
    ---------
    1. Retrieve density `rho` and the element coordinates.
    2. Compute the element volume from `det(J) / 6`.
    3. Generate the analytical consistent mass scaling if not lumped.
    4. Distribute mass diagonally if `lumped=True`.
    """
    Xe = as_float_array(Xe)
    props = as_float_array(Ge).reshape(-1)
    rho = props[2] if props.size > 2 else 1.0

    # Volume from Jacobian
    dN = np.array(
        [[1.0, 0.0, 0.0, -1.0], [0.0, 1.0, 0.0, -1.0], [0.0, 0.0, 1.0, -1.0]],
        dtype=float,
    )
    J = dN @ Xe
    V = abs(np.linalg.det(J)) / 6.0

    if lumped:
        return (rho * V / 4.0) * np.eye(12, dtype=float)

    I3 = np.eye(3, dtype=float)
    scalar = np.array(
        [[2, 1, 1, 1], [1, 2, 1, 1], [1, 1, 2, 1], [1, 1, 1, 2]],
        dtype=float,
    )
    Me = (rho * V / 20.0) * np.kron(scalar, I3)
    return Me


def mT4e(M, T, X, G, *, lumped: bool = False):
    """
    Assemble T4 element mass matrices into the global mass matrix.

    Parameters
    ----------
    M : ndarray or sparse, shape (ndof, ndof)
        Global mass matrix (modified in place).
    T : array_like, shape (nel, 5)
        Topology ``[n1, n2, n3, n4, mat_id]``.
    X : array_like, shape (nn, 3)
        Nodal coordinates.
    G : array_like
        Material table ``[E, nu, rho]``.
    lumped : bool
        If True, assemble lumped mass.

    Returns
    -------
    M : ndarray or sparse
        Updated global mass matrix.

    Mathematical Formulation
    ------------------------
    Global mass matrix assembly combines the T4 element mass matrices in the
    standard form `L_e.T @ Me @ L_e`.

    Algorithm
    ---------
    1. Parse global nodal coordinates and topology mapping.
    2. Extract density, compute batched volumes via determinant.
    3. Generate consistent or lumped elemental mass matrices.
    4. Scatter the element mass properties into the global system matrix.
    """
    topology = as_float_array(T)
    coordinates = as_float_array(X)
    nodes = topology[:, :4].astype(int) - 1
    materials = as_float_array(G)[topology[:, -1].astype(int) - 1]

    Xe = coordinates[nodes]
    dN_ref = np.array(
        [[1.0, 0.0, 0.0, -1.0], [0.0, 1.0, 0.0, -1.0], [0.0, 0.0, 1.0, -1.0]],
        dtype=float,
    )
    J = np.einsum("ij,ejk->eik", dN_ref, Xe)
    volumes = np.abs(np.linalg.det(J)) / 6.0
    rho = materials[:, 2] if materials.shape[1] > 2 else np.ones(topology.shape[0])

    indices = element_dof_indices(nodes, 3, one_based=False)

    if lumped:
        mass_per_node = rho * volumes / 4.0
        for e in range(topology.shape[0]):
            idx = indices[e]
            for k in idx:
                M[k, k] += mass_per_node[e]
        return M

    I3 = np.eye(3, dtype=float)
    scalar = np.array(
        [[2, 1, 1, 1], [1, 2, 1, 1], [1, 1, 2, 1], [1, 1, 1, 2]],
        dtype=float,
    )
    block = np.kron(scalar, I3)
    factors = rho * volumes / 20.0
    element_matrices = factors[:, None, None] * block[None, :, :]

    if is_sparse(M) and sp is not None:
        scatter_rows = np.broadcast_to(
            indices[:, :, None], element_matrices.shape
        ).reshape(-1)
        scatter_cols = np.broadcast_to(
            indices[:, None, :], element_matrices.shape
        ).reshape(-1)
        delta = sp.coo_matrix(
            (element_matrices.reshape(-1), (scatter_rows, scatter_cols)),
            shape=M.shape,
            dtype=float,
        )
        return (M.tocsr() + delta.tocsr()).tolil()
    np.add.at(M, (indices[:, :, None], indices[:, None, :]), element_matrices)
    return M


def meh8e(Xe, Ge, *, lumped: bool = False):
    """
    Compute the element mass matrix for an 8-node hexahedral solid element.

    Consistent (24x24) via 2x2x2 Gauss quadrature:
        M = sum_gp  rho * N^T N * det(J) * w

    Lumped via HRZ (diagonal scaling to preserve total mass).

    Parameters
    ----------
    Xe : array_like, shape (8, 3)
        Nodal coordinates.
    Ge : array_like
        Material row ``[E, nu, rho]``.
    lumped : bool
        If True, return diagonally lumped mass.

    Returns
    -------
    Me : ndarray, shape (24, 24)

    Mathematical Formulation
    ------------------------
    Using standard `2 x 2 x 2` Gauss integration, the consistent mass matrix
    is obtained by summing `rho * N_i.T @ N_i * det(J_i)` over the eight
    Gauss points. The lumped version uses HRZ-style diagonal scaling to
    preserve the total element mass.

    Algorithm
    ---------
    1. Initialize the elemental mass accumulator tensor.
    2. Integrate density over 8 Gaussian nodes using shape functions.
    3. Scale density with the Jacobian determinant per volume step.
    4. If lumped, apply HRZ extraction along the diagonal elements.
    """
    Xe = as_float_array(Xe)
    props = as_float_array(Ge).reshape(-1)
    rho = props[2] if props.size > 2 else 1.0

    gauss_points = np.array(
        [
            [-1, -1, -1],
            [1, -1, -1],
            [1, 1, -1],
            [-1, 1, -1],
            [-1, -1, 1],
            [1, -1, 1],
            [1, 1, 1],
            [-1, 1, 1],
        ],
        dtype=float,
    ) / np.sqrt(3.0)

    Me = np.zeros((24, 24), dtype=float)

    for gp in gauss_points:
        ri, rj, rk = gp
        # Shape functions for H8
        N_vals = (
            np.array(
                [
                    (1 - ri) * (1 - rj) * (1 - rk),
                    (1 + ri) * (1 - rj) * (1 - rk),
                    (1 + ri) * (1 + rj) * (1 - rk),
                    (1 - ri) * (1 + rj) * (1 - rk),
                    (1 - ri) * (1 - rj) * (1 + rk),
                    (1 + ri) * (1 - rj) * (1 + rk),
                    (1 + ri) * (1 + rj) * (1 + rk),
                    (1 - ri) * (1 + rj) * (1 + rk),
                ],
                dtype=float,
            )
            / 8.0
        )

        N_mat = np.zeros((3, 24), dtype=float)
        N_mat[0, 0::3] = N_vals
        N_mat[1, 1::3] = N_vals
        N_mat[2, 2::3] = N_vals

        dN = _hexa_dN(ri, rj, rk)
        Jt = dN @ Xe
        detJ = np.linalg.det(Jt)
        Me += rho * (N_mat.T @ N_mat) * detJ  # weights are 1 for 2x2x2

    if lumped:
        diag_vals = np.diag(Me).copy()
        total_mass = Me.sum() / 3.0  # total mass = trace / ndim
        diag_sum = diag_vals.sum()
        if diag_sum > 0:
            diag_vals *= (total_mass * 3.0) / diag_sum
        return np.diag(diag_vals)

    return Me


def mh8e(M, T, X, G, *, lumped: bool = False):
    """
    Assemble H8 element mass matrices into the global mass matrix.

    Parameters
    ----------
    M : ndarray or sparse, shape (ndof, ndof)
        Global mass matrix (modified in place).
    T : array_like, shape (nel, 9)
        Topology ``[n1, ..., n8, mat_id]``.
    X : array_like, shape (nn, 3)
        Nodal coordinates.
    G : array_like
        Material table ``[E, nu, rho]``.
    lumped : bool
        If True, assemble lumped mass.

    Returns
    -------
    M : ndarray or sparse
        Updated global mass matrix.

    Mathematical Formulation
    ------------------------
    Global mass assembly adds each element matrix in the standard form
    `L_e.T @ Me @ L_e`.

    Algorithm
    ---------
    1. Query global coordinates and corresponding element connectivity.
    2. Iteratively process individual element structural arrays.
    3. Determine the corresponding `24 x 24` element mass entries.
    4. Combine and scatter locally computed coefficients into the global mass.
    """
    topology = as_float_array(T)
    coordinates = as_float_array(X)
    nodes = topology[:, :8].astype(int) - 1
    materials = as_float_array(G)[topology[:, -1].astype(int) - 1]

    indices = element_dof_indices(nodes, 3, one_based=False)

    for e in range(topology.shape[0]):
        Xe = coordinates[nodes[e]]
        Me = meh8e(Xe, materials[e], lumped=lumped)
        idx = indices[e]
        if is_sparse(M):
            M[np.ix_(idx, idx)] = M[np.ix_(idx, idx)] + Me
        else:
            M[np.ix_(idx, idx)] += Me
    return M


__all__ = [
    "keT4e",
    "keh8e",
    "kT4e",
    "kh8e",
    "meT4e",
    "meh8e",
    "mT4e",
    "mh8e",
    "qeT4e",
    "qeh8e",
    "qT4e",
    "qh8e",
]
