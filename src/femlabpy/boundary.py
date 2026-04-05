"""
Boundary-condition application and constrained linear-system solvers.

Workflow role
-------------
This module sits between assembly and solution. ``setbc`` applies prescribed
displacements by direct elimination, while ``solve_lag_general`` and
``solve_lag`` build augmented systems for exact linear constraints using
Lagrange multipliers.

What to read first
------------------
- ``setbc`` for the standard fixed-support workflow.
- ``solve_lag_general`` for constraint equations such as periodicity or multi-
  point conditions.
- ``solve_lag`` for the legacy wrapper that keeps the original FemLab calling
  convention.
"""

from __future__ import annotations

import numpy as np

from ._helpers import (
    as_float_array,
    is_sparse,
    max_abs_diagonal,
    solve_linear_system,
)

try:
    import scipy.sparse as sp
except ImportError:  # pragma: no cover
    sp = None


def solve_lag_general(
    K,
    p,
    G,
    Q=None,
    *,
    scale: float | None = None,
    return_lagrange: bool = False,
):
    """Solve a linear system with general linear constraints ``G u = Q``.

    The augmented system is scaled to keep the constraint rows numerically
    compatible with the stiffness matrix, matching the legacy toolbox pattern.

    Algorithm
    ---------
    1. Define the scaling factor $\alpha = 10^{-2} \times \max_i |K_{ii}|$.
    2. Form the scaled constraint matrix $\bar{\mathbf{G}} = \alpha \mathbf{G}$ and right-hand side $\bar{\mathbf{Q}} = \alpha \mathbf{Q}$.
    3. Assemble the augmented matrix:
       $$
       \bar{\mathbf{K}} = \begin{bmatrix} \mathbf{K} & \bar{\mathbf{G}}^T \\ \bar{\mathbf{G}} & \mathbf{0} \end{bmatrix}
       $$
    4. Assemble the augmented right-hand side:
       $$
       \bar{\mathbf{p}} = \begin{bmatrix} \mathbf{p} \\ \bar{\mathbf{Q}} \end{bmatrix}
       $$
    5. Solve the augmented system $\bar{\mathbf{K}} \mathbf{x} = \bar{\mathbf{p}}$ to obtain the solution $\mathbf{u}$ and scaled Lagrange multipliers $\bar{\boldsymbol{\lambda}}$.
    6. If requested, recover the physical multipliers $\boldsymbol{\lambda} = \alpha \bar{\boldsymbol{\lambda}}$.
    """

    constraint_matrix = as_float_array(G)
    if constraint_matrix.size == 0:
        solution = solve_linear_system(K, p)
        if return_lagrange:
            return solution, np.zeros((0, 1), dtype=float)
        return solution

    if constraint_matrix.ndim == 1:
        constraint_matrix = constraint_matrix.reshape(1, -1)

    system_size = K.shape[0]
    if constraint_matrix.shape[1] != system_size:
        raise ValueError(
            "Constraint matrix width must match the number of system DOFs."
        )

    if Q is None:
        constraint_rhs = np.zeros((constraint_matrix.shape[0], 1), dtype=float)
    else:
        constraint_rhs = as_float_array(Q).reshape(-1, 1)
        if constraint_rhs.shape[0] != constraint_matrix.shape[0]:
            raise ValueError(
                "Constraint RHS height must match the number of constraint rows."
            )

    if scale is None:
        scale = 1.0e-2 * max_abs_diagonal(K)
        if scale == 0.0:
            scale = 1.0

    Gbar = scale * constraint_matrix
    Qbar = scale * constraint_rhs

    if is_sparse(K) and sp is not None:
        Kbar = sp.bmat(
            [
                [K, sp.csr_matrix(Gbar.T)],
                [sp.csr_matrix(Gbar), None],
            ],
            format="csr",
        )
    else:
        Kbar = np.block(
            [
                [as_float_array(K), Gbar.T],
                [
                    Gbar,
                    np.zeros(
                        (constraint_matrix.shape[0], constraint_matrix.shape[0]),
                        dtype=float,
                    ),
                ],
            ]
        )

    pbar = np.vstack([as_float_array(p).reshape(-1, 1), Qbar])
    augmented = solve_linear_system(Kbar, pbar)
    solution = augmented[:system_size]

    if return_lagrange:
        lagrange = augmented[system_size:] * scale
        return solution, lagrange
    return solution


def setbc(K, p, C, dof: int = 1):
    """
    Apply boundary conditions using direct elimination.

    For each constrained DOF *j* with prescribed value *d*:

    1. Transfer the coupling forces to the RHS:  ``p -= K[:, j] * d``
    2. Zero out the row and column:  ``K[j, :] = 0``, ``K[:, j] = 0``
    3. Place a spring stiffness on the diagonal:  ``K[j, j] = ks``
    4. Set the load entry:  ``p[j] = ks * d``

    This is the standard textbook direct-elimination approach matching
    the Scilab FemLab ``setbc.sci`` convention (row/column zeroing with
    ``ks = 0.1 * max(diag(K))``), with the additional correction that
    coupling forces are transferred to the RHS *before* zeroing, so
    non-zero prescribed displacements are handled correctly.

    Algorithm
    ---------
    1. Calculate a penalty stiffness $k_s = 0.1 \times \max_i |K_{ii}|$.
    2. For each constraint $j$ with prescribed value $d$:
       a. Transfer coupling forces to the RHS: $p_i \leftarrow p_i - K_{ij} d$ for $i \neq j$.
       b. Zero out the corresponding row and column: $K_{ij} = 0$ and $K_{ji} = 0$ for all $i$.
       c. Place the penalty stiffness on the diagonal: $K_{jj} = k_s$.
       d. Set the corresponding load entry: $p_j = k_s d$.

    Parameters
    ----------
    K : ndarray or sparse matrix, shape (ndof, ndof)
        Global stiffness matrix (modified in place).

    p : ndarray, shape (ndof, 1)
        Load vector (modified in place).

    C : array_like, shape (nbc, 2) or (nbc, 3)
        Boundary condition array:
        - For dof=1: each row is [node, value]
        - For dof>1: each row is [node, local_dof, value]
        Node indices are 1-based.

    dof : int, default=1
        Number of DOFs per node.
        Use dof=2 for 2D problems, dof=3 for 3D problems.

    Returns
    -------
    K : ndarray or sparse matrix
        Modified stiffness matrix.

    p : ndarray
        Modified load vector.

    ks : float
        Spring stiffness used for the constrained DOFs.

    Examples
    --------
    >>> from femlabpy import init, setbc
    >>> K, p, q = init(nn=10, dof=2)
    >>> # Fix node 1 (both DOFs) and node 2 (x-direction only)
    >>> C = np.array([
    ...     [1, 1, 0.0],  # node 1, ux = 0
    ...     [1, 2, 0.0],  # node 1, uy = 0
    ...     [2, 1, 0.0],  # node 2, ux = 0
    ... ])
    >>> K, p, _ = setbc(K, p, C, dof=2)
    """
    constraints = as_float_array(C)
    if constraints.size == 0:
        return K, p, 0.0

    ks = 0.1 * max_abs_diagonal(K)
    if ks == 0.0:
        ks = 1.0

    if dof == 1:
        cdofs = constraints[:, 0].astype(int) - 1
        cvals = constraints[:, -1]
    else:
        cdofs = (
            (constraints[:, 0].astype(int) - 1) * dof
            + constraints[:, 1].astype(int)
            - 1
        )
        cvals = constraints[:, -1]

    # Convert sparse matrices to lil for efficient row/column zeroing.
    sparse = is_sparse(K)
    if sparse:
        K = K.tolil()

    for k in range(len(cdofs)):
        j = int(cdofs[k])
        val = cvals[k]
        # Transfer coupling forces to RHS *before* zeroing the column.
        if val != 0.0:
            if sparse:
                col_j = np.asarray(K[:, j].toarray()).ravel()
            else:
                col_j = K[:, j].copy()
            p[:, 0] -= col_j * val
        # Zero row and column, set diagonal spring.
        K[j, :] = 0
        K[:, j] = 0
        K[j, j] = ks
        p[j, 0] = ks * val

    return K, p, float(ks)


def solve_lag(K, p, C=None, dof: int = 1, *, return_lagrange: bool = False):
    """
    Solve a linear system with nodal Dirichlet constraints via Lagrange multipliers.

    Parameters
    ----------
    K, p:
        Linear system ``K u = p`` before applying displacement constraints.
    C:
        Legacy constraint table. For scalar problems, each row is
        ``[node, value]``. For vector problems, each row is
        ``[node, local_dof, value]``.
    dof:
        Degrees of freedom per node.
    return_lagrange:
        When ``True``, also return the recovered Lagrange multipliers.

    Algorithm
    ---------
    1. Map the nodal constraints given in $\mathbf{C}$ to global degree of freedom indices $j$.
    2. Construct the Boolean constraint matrix $\mathbf{G}$ where $G_{ij} = 1$ if constraint $i$ acts on DOF $j$, and $0$ otherwise.
    3. Extract the constraint values into the vector $\mathbf{Q}$.
    4. Delegate the solution to the general augmented Lagrangian solver to compute $\mathbf{u}$ and $\boldsymbol{\lambda}$ satisfying:
       $$
       \begin{bmatrix} \mathbf{K} & \mathbf{G}^T \\ \mathbf{G} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \mathbf{u} \\ \boldsymbol{\lambda} \end{bmatrix} = \begin{bmatrix} \mathbf{p} \\ \mathbf{Q} \end{bmatrix}
       $$

    Returns
    -------
    ndarray or tuple[ndarray, ndarray]
        Constrained solution vector, optionally paired with the multiplier
        vector.
    """
    if C is None:
        solution = solve_linear_system(K, p)
        if return_lagrange:
            return solution, np.zeros((0, 1), dtype=float)
        return solution
    constraints = as_float_array(C)
    if constraints.size == 0:
        solution = solve_linear_system(K, p)
        if return_lagrange:
            return solution, np.zeros((0, 1), dtype=float)
        return solution

    n_constraints = constraints.shape[0]
    system_size = K.shape[0]

    G = np.zeros((n_constraints, system_size), dtype=float)
    Q = constraints[:, -1].reshape(-1, 1)
    if dof == 1:
        indices = constraints[:, 0].astype(int) - 1
    else:
        indices = (
            (constraints[:, 0].astype(int) - 1) * dof
            + constraints[:, 1].astype(int)
            - 1
        )
    G[np.arange(n_constraints), indices] = 1.0

    return solve_lag_general(K, p, G, Q, return_lagrange=return_lagrange)


def rnorm(f, C, dof: int):
    """
    Return the residual norm restricted to unconstrained degrees of freedom.

    Parameters
    ----------
    f:
        Residual or force vector.
    C:
        Legacy constraint table ``[node, local_dof, value]``.
    dof:
        Degrees of freedom per node.

    Algorithm
    ---------
    1. Map the constrained nodal degrees of freedom in $\mathbf{C}$ to global indices $I_c$.
    2. Identify the unconstrained index set $I_u = \{1, 2, \ldots, N\} \setminus I_c$.
    3. Compute the Euclidean norm over the unconstrained components:
       $$
       \|\mathbf{f}_{I_u}\|_2 = \sqrt{\sum_{i \in I_u} f_i^2}
       $$

    Returns
    -------
    float
        Euclidean norm of the unconstrained residual entries.
    """
    force = as_float_array(f).reshape(-1)
    constraints = as_float_array(C)
    fixed = np.zeros(force.shape[0], dtype=bool)
    indices = (
        (constraints[:, 0].astype(int) - 1) * dof + constraints[:, 1].astype(int) - 1
    )
    fixed[indices] = True
    return float(np.linalg.norm(force[~fixed]))


__all__ = ["rnorm", "setbc", "solve_lag", "solve_lag_general"]
