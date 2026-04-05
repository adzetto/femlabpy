"""
Global array allocation and small indexing helpers for the femlabpy workflow.

Workflow role
-------------
This module is the normal starting point after the node table and the number of
degrees of freedom per node are known. It creates the zeroed global arrays that
the assembly, boundary, modal, and dynamics modules fill later in the solve
sequence.

Public entry points
-------------------
- ``init`` allocates stiffness, optional mass, external load, and internal
  force arrays with dense or sparse storage.
- ``rows`` and ``cols`` are lightweight shape helpers re-exported from the
  private helper layer because the legacy translation code uses them often.
"""

from __future__ import annotations

from ._helpers import cols, rows, zeros_matrix, zeros_vector


def init(nn: int, dof: int, *, dynamic: bool = False, use_sparse: bool | None = None):
    """
    Initialize FEM arrays for a problem with nn nodes and dof DOFs per node.

    Parameters
    ----------
    nn : int
        Number of nodes in the mesh.

    dof : int
        Degrees of freedom per node (2 for 2D, 3 for 3D).

    dynamic : bool, default False
        If True, also allocate a mass matrix ``M`` and return it.

    use_sparse : bool, optional
        If True, use scipy.sparse.lil_matrix for K (and M).
        If False, use dense numpy array.
        If None (default), automatically use sparse for nn >= 1000.

    Returns
    -------
    K : ndarray or sparse matrix, shape (ndof, ndof)
        Global stiffness matrix (initialized to zero).

    M : ndarray or sparse matrix, shape (ndof, ndof)
        Global mass matrix (only when ``dynamic=True``).

    p : ndarray, shape (ndof, 1)
        Load vector (initialized to zero).

    q : ndarray, shape (ndof, 1)
        Internal force vector (initialized to zero).

    Notes
    -----
    Total DOFs = nn * dof

    Algorithm
    ---------
    1. Compute the total degrees of freedom $N = \text{nn} \times \text{dof}$.
    2. Allocate a zero matrix $\mathbf{K} \in \mathbb{R}^{N \times N}$ for the stiffness matrix.
    3. Allocate zero vectors $\mathbf{p}, \mathbf{q} \in \mathbb{R}^{N \times 1}$ for the load and internal force vectors.
    4. If $\text{dynamic}$ is `True`, also allocate a mass matrix $\mathbf{M} \in \mathbb{R}^{N \times N}$.

    When ``dynamic=False`` (default), returns ``(K, p, q)`` for backward
    compatibility.  When ``dynamic=True``, returns ``(K, M, p, q)``.

    Examples
    --------
    >>> from femlabpy import init
    >>> K, p, q = init(nn=100, dof=2)  # 100 nodes, 2D problem
    >>> K.shape
    (200, 200)

    >>> # Dynamic allocation with mass matrix
    >>> K, M, p, q = init(nn=100, dof=2, dynamic=True)

    >>> # Force sparse storage
    >>> K, p, q = init(nn=50, dof=2, use_sparse=True)
    """
    total_dofs = int(nn) * int(dof)
    if use_sparse is None:
        use_sparse = nn >= 1000
    stiffness = zeros_matrix(total_dofs, use_sparse=use_sparse)
    load = zeros_vector(total_dofs)
    internal = zeros_vector(total_dofs)
    if dynamic:
        mass = zeros_matrix(total_dofs, use_sparse=use_sparse)
        return stiffness, mass, load, internal
    return stiffness, load, internal


__all__ = ["cols", "init", "rows"]
