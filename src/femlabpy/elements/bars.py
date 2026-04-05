"""
Geometrically nonlinear bar and truss element kernels.

Workflow role
-------------
This module contains both element-level routines and assembled global wrappers
for bar-like members. It is the main reference for one-dimensional axial
behavior, geometric nonlinearity, and the simplest mass-matrix construction in
the codebase.

Public entry points
-------------------
- ``kebar`` and ``qebar`` work at single-element level.
- ``kbar`` and ``qbar`` assemble those contributions over a full topology table.
- ``mebar`` and ``mbar`` provide consistent or lumped mass matrices for
  dynamic analyses.
"""

from __future__ import annotations

import numpy as np

from .._helpers import (
    as_float_array,
    cols,
    element_dof_indices,
    is_sparse,
    rows,
)

try:
    import scipy.sparse as sp
except ImportError:  # pragma: no cover
    sp = None


def kebar(Xe0, Xe1, Ge):
    """Compute the tangent stiffness matrix of a geometrically nonlinear bar element.

    Mathematical Formulation
    ------------------------
    The tangent stiffness matrix $K$ is the sum of the material stiffness matrix $K_m$ and the geometric stiffness matrix $K_g$:
    $$K = K_m + K_g$$
    where
    $$K_m = \\frac{E A}{L_0^3} \\begin{bmatrix} a_1 a_1^T & -a_1 a_1^T \\\\ -a_1 a_1^T & a_1 a_1^T \\end{bmatrix}$$
    $$K_g = \\frac{N}{L_0} \\begin{bmatrix} I & -I \\\\ -I & I \\end{bmatrix}$$
    Here, $a_1$ is the current element vector, $L_0$ is the initial length, $A$ is the cross-sectional area, $E$ is the Young's modulus, and $N$ is the axial force.

    Algorithm
    ---------
    1. Compute the initial length $L_0$ and current length $L_1$.
    2. Compute the Green-Lagrange strain $\\epsilon = \\frac{L_1^2 - L_0^2}{2 L_0^2}$.
    3. Compute the normal force $N = E A \\epsilon$.
    4. Assemble and return the tangent stiffness matrix $K = K_m + K_g$.
    """
    initial = as_float_array(Xe0)
    current = as_float_array(Xe1)
    props = as_float_array(Ge).reshape(-1)

    a0 = (initial[1] - initial[0]).reshape(-1, 1)
    l0 = float(np.linalg.norm(a0))
    a1 = (current[1] - current[0]).reshape(-1, 1)
    l1 = float(np.linalg.norm(a1))

    A = props[0]
    E = props[1] if props.size > 1 else 1.0
    strain = 0.5 * (l1**2 - l0**2) / l0**2
    normal_force = A * E * strain
    identity = np.eye(a0.shape[0], dtype=float)
    return (E * A / l0**3) * np.block(
        [[a1 @ a1.T, -a1 @ a1.T], [-a1 @ a1.T, a1 @ a1.T]]
    ) + (normal_force / l0) * np.block([[identity, -identity], [-identity, identity]])


def qebar(Xe0, Xe1, Ge):
    """Compute the internal-force response of a single geometrically nonlinear bar.

    Mathematical Formulation
    ------------------------
    The internal force vector $q$ for a bar element is derived from the variation of the strain energy. It can be expressed as:
    $$q = \\frac{A \\sigma}{L_0} \\begin{bmatrix} -a_1 \\\\ a_1 \\end{bmatrix}$$
    where $A$ is the cross-sectional area, $\\sigma$ is the normal stress, $L_0$ is the initial length, and $a_1$ is the current element vector.
    The Green-Lagrange strain is given by $\\epsilon = \\frac{L_1^2 - L_0^2}{2 L_0^2}$ and the stress is $\\sigma = E \\epsilon$.

    Algorithm
    ---------
    1. Compute the initial and current element vectors and lengths ($L_0$, $L_1$).
    2. Compute the Green-Lagrange strain $\\epsilon$.
    3. Compute the normal stress $\\sigma = E \\epsilon$.
    4. Calculate the internal force vector $q$.
    5. Return the internal force vector $q$, stress $\\sigma$, and strain $\\epsilon$.
    """
    initial = as_float_array(Xe0)
    current = as_float_array(Xe1)
    props = as_float_array(Ge).reshape(-1)

    a0 = (initial[1] - initial[0]).reshape(-1, 1)
    l0 = float(np.linalg.norm(a0))
    a1 = (current[1] - current[0]).reshape(-1, 1)
    l1 = float(np.linalg.norm(a1))

    A = props[0]
    E = props[1] if props.size > 1 else 1.0
    strain = 0.5 * (l1**2 - l0**2) / l0**2
    stress = E * strain
    qe = (A * stress / l0) * np.vstack([-a1, a1])
    return qe, float(stress), float(strain)


def kbar(K, T, X, G, u=None):
    """Assemble bar or truss tangent stiffness contributions into the global matrix.

    Mathematical Formulation
    ------------------------
    The global tangent stiffness matrix $\\mathbf{K}$ is updated by assembling the contributions of each geometrically nonlinear bar element. For each element $e$, the stiffness matrix $K^{(e)}$ is:
    $$K^{(e)} = \\frac{E A}{L_0^3} \\begin{bmatrix} a_1 a_1^T & -a_1 a_1^T \\\\ -a_1 a_1^T & a_1 a_1^T \\end{bmatrix} + \\frac{N}{L_0} \\begin{bmatrix} I & -I \\\\ -I & I \\end{bmatrix}$$
    The global stiffness matrix is obtained by summing the element stiffness matrices using the connectivity matrix $\\mathbf{A}^{(e)}$:
    $$\\mathbf{K} = \\sum_{e=1}^{N_{el}} (\\mathbf{A}^{(e)})^T K^{(e)} \\mathbf{A}^{(e)}$$

    Algorithm
    ---------
    1. Determine the current nodal coordinates based on initial coordinates $X$ and displacements $u$.
    2. Compute the initial and current lengths of all bar elements.
    3. Calculate the strains and normal forces for all elements.
    4. Construct the tangent stiffness matrices (material and geometric) for all elements in a vectorized manner.
    5. Assemble the element matrices into the global stiffness matrix $K$ using sparse or dense operations.
    6. Return the updated global matrix $K$.
    """
    X = as_float_array(X)
    topology = as_float_array(T)
    if u is None:
        current = X
    else:
        current = X + as_float_array(u).reshape(rows(X), cols(X))
    element_nodes = topology[:, :-1].astype(int) - 1
    props = as_float_array(G)[topology[:, -1].astype(int) - 1]
    initial = X[element_nodes]
    current_nodes = current[element_nodes]
    a0 = initial[:, 1, :] - initial[:, 0, :]
    a1 = current_nodes[:, 1, :] - current_nodes[:, 0, :]
    l0 = np.linalg.norm(a0, axis=1)
    l1 = np.linalg.norm(a1, axis=1)
    area = props[:, 0]
    modulus = props[:, 1] if props.shape[1] > 1 else np.ones_like(area)
    strain = 0.5 * (l1**2 - l0**2) / l0**2
    normal_force = area * modulus * strain
    a1a1 = np.einsum("ei,ej->eij", a1, a1)
    identity = np.eye(cols(X), dtype=float)[None, :, :]
    axial = (modulus * area / l0**3)[:, None, None] * a1a1
    geometric = (normal_force / l0)[:, None, None] * identity
    upper = np.concatenate([axial + geometric, -(axial + geometric)], axis=2)
    lower = np.concatenate([-(axial + geometric), axial + geometric], axis=2)
    element_matrices = np.concatenate([upper, lower], axis=1)
    indices = element_dof_indices(element_nodes, cols(X), one_based=False)
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


def qbar(q, T, X, G, u=None):
    """Assemble bar or truss internal forces and element output quantities.

    Mathematical Formulation
    ------------------------
    The global internal force vector $\\mathbf{q}$ is updated by assembling the local internal force vectors $q^{(e)}$ of all elements. For each element $e$:
    $$q^{(e)} = \\frac{A \\sigma}{L_0} \\begin{bmatrix} -a_1 \\\\ a_1 \\end{bmatrix}$$
    The global internal forces are computed by:
    $$\\mathbf{q} = \\sum_{e=1}^{N_{el}} (\\mathbf{A}^{(e)})^T q^{(e)}$$

    Algorithm
    ---------
    1. Evaluate current coordinates from initial coordinates and displacements.
    2. Compute initial lengths, current lengths, strains, and stresses for all elements.
    3. Compute local internal force vectors for all elements.
    4. Assemble the local internal force vectors into the global internal force vector $q$.
    5. Return the updated global internal force vector $q$, and the arrays of stresses and strains.
    """
    X = as_float_array(X)
    topology = as_float_array(T)
    if u is None:
        current = X
    else:
        current = X + as_float_array(u).reshape(rows(X), cols(X))
    element_nodes = topology[:, :-1].astype(int) - 1
    props = as_float_array(G)[topology[:, -1].astype(int) - 1]
    initial = X[element_nodes]
    current_nodes = current[element_nodes]
    a0 = initial[:, 1, :] - initial[:, 0, :]
    a1 = current_nodes[:, 1, :] - current_nodes[:, 0, :]
    l0 = np.linalg.norm(a0, axis=1)
    l1 = np.linalg.norm(a1, axis=1)
    area = props[:, 0]
    modulus = props[:, 1] if props.shape[1] > 1 else np.ones_like(area)
    strain = (0.5 * (l1**2 - l0**2) / l0**2).reshape(-1, 1)
    stress = (modulus[:, None] * strain).reshape(-1, 1)
    element_vectors = (area * stress[:, 0] / l0)[:, None] * np.concatenate(
        [-a1, a1], axis=1
    )
    indices = element_dof_indices(element_nodes, cols(X), one_based=False)
    np.add.at(q[:, 0], indices.reshape(-1), element_vectors.reshape(-1))
    return q, stress, strain


def mebar(Xe, Ge, dof: int = 2, *, lumped: bool = False):
    """
    Compute the element mass matrix for a 2-node bar/truss element.

    Consistent mass:
        M = (rho * A * L / 6) * [[2, 1], [1, 2]] tensor I_dof

    Lumped mass:
        M = (rho * A * L / 2) * I_{2*dof}

    Mathematical Formulation
    ------------------------
    The consistent mass matrix $M$ is derived using the shape functions $N_i$:
    $$M = \\int_0^L \\rho A N^T N \\, dx = \\frac{\\rho A L}{6} \\begin{bmatrix} 2 & 1 \\\\ 1 & 2 \\end{bmatrix} \\otimes I_{dof}$$
    The lumped mass matrix employs a diagonal representation, dividing the total mass $\\rho A L$ equally among the nodes:
    $$M_{lumped} = \\frac{\\rho A L}{2} I_{2 \\cdot dof}$$

    Algorithm
    ---------
    1. Extract the area $A$, density $\\rho$, and initial coordinates to find length $L$.
    2. If a lumped mass matrix is requested, compute and return the diagonal matrix.
    3. If a consistent mass matrix is requested, construct the block matrix using the standard formulation.
    4. Return the element mass matrix $M$.

    Parameters
    ----------
    Xe : array_like, shape (2, ndim)
        Initial (undeformed) nodal coordinates.
    Ge : array_like
        Material row ``[A, E, rho]``.  If ``rho`` is omitted it defaults to 1.
    dof : int
        Degrees of freedom per node.
    lumped : bool
        If True, return the diagonally lumped mass matrix.

    Returns
    -------
    Me : ndarray, shape (2*dof, 2*dof)
    """
    Xe = as_float_array(Xe)
    props = as_float_array(Ge).reshape(-1)
    A = props[0]
    rho = props[2] if props.size > 2 else 1.0
    a0 = Xe[1] - Xe[0]
    L = float(np.linalg.norm(a0))
    size = 2 * dof

    if lumped:
        return (rho * A * L / 2.0) * np.eye(size, dtype=float)

    # Consistent: (rho*A*L/6) * [[2I, 1I],[1I, 2I]]
    I_d = np.eye(dof, dtype=float)
    Me = (rho * A * L / 6.0) * np.block(
        [[2.0 * I_d, 1.0 * I_d], [1.0 * I_d, 2.0 * I_d]]
    )
    return Me


def mbar(M, T, X, G, dof: int = 2, *, lumped: bool = False):
    """
    Assemble bar/truss mass matrices into the global mass matrix.

    Mathematical Formulation
    ------------------------
    The global mass matrix $\\mathbf{M}$ is constructed by assembling the local mass matrices $M^{(e)}$ of all elements.
    $$\\mathbf{M} = \\sum_{e=1}^{N_{el}} (\\mathbf{A}^{(e)})^T M^{(e)} \\mathbf{A}^{(e)}$$
    Depending on the formulation, $M^{(e)}$ is either the consistent or lumped mass matrix.

    Algorithm
    ---------
    1. Extract initial nodal coordinates and properties for all elements.
    2. Compute the length $L$, area $A$, and density $\\rho$ for each element.
    3. If lumped mass is specified, add half the total element mass to the diagonal entries of the corresponding global degrees of freedom.
    4. If consistent mass is specified, construct the full element mass matrices and assemble them into the global mass matrix $M$.
    5. Return the updated global mass matrix $M$.

    Parameters
    ----------
    M : ndarray or sparse, shape (ndof, ndof)
        Global mass matrix (modified in place).
    T : array_like, shape (nel, 3)
        Topology table ``[n1, n2, mat_id]``.
    X : array_like, shape (nn, ndim)
        Nodal coordinates.
    G : array_like
        Material table with rows ``[A, E, rho]``.
    dof : int
        Degrees of freedom per node.
    lumped : bool
        If True, assemble lumped mass.

    Returns
    -------
    M : ndarray or sparse
        Updated global mass matrix.
    """
    X = as_float_array(X)
    topology = as_float_array(T)
    element_nodes = topology[:, :-1].astype(int) - 1
    props = as_float_array(G)[topology[:, -1].astype(int) - 1]

    initial = X[element_nodes]
    a0 = initial[:, 1, :] - initial[:, 0, :]
    lengths = np.linalg.norm(a0, axis=1)
    area = props[:, 0]
    rho = props[:, 2] if props.shape[1] > 2 else np.ones_like(area)

    indices = element_dof_indices(element_nodes, dof, one_based=False)

    if lumped:
        mass_per_node = 0.5 * rho * area * lengths
        for e in range(topology.shape[0]):
            idx = indices[e]
            for k in idx:
                M[k, k] += mass_per_node[e]
    else:
        I_d = np.eye(dof, dtype=float)
        block_22 = np.block([[2.0 * I_d, 1.0 * I_d], [1.0 * I_d, 2.0 * I_d]])
        factors = rho * area * lengths / 6.0
        element_matrices = factors[:, None, None] * block_22[None, :, :]
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


__all__ = ["kbar", "kebar", "mbar", "mebar", "qbar", "qebar"]
