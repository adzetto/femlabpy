"""
Scatter-assembly helpers that map element-level results into global arrays.

Workflow role
-------------
Element kernels in :mod:`femlabpy.elements` compute local matrices and vectors
for one element at a time. This module converts the one-based topology row into
global degree-of-freedom indices and adds those local contributions into the
global stiffness matrix or internal-force vector in place.

Reading guide
-------------
Read ``assmk`` first if you want to understand matrix assembly and ``assmq``
second for force-vector assembly. Both functions share the same indexing path,
so once the topology-to-DOF mapping is clear the rest of the global workflow is
much easier to follow.
"""

from __future__ import annotations

import numpy as np

from ._helpers import as_float_array, is_sparse, node_dof_indices, topology_nodes


def assmk(K, Ke, Te, dof: int = 1):
    """
    Assemble one element stiffness matrix into the global stiffness matrix.

    Parameters
    ----------
    K:
        Global stiffness matrix updated in place.
    Ke:
        Element stiffness matrix ordered by the element's local DOFs.
    Te:
        Legacy FemLab topology row whose last entry stores the property id.
    dof:
        Degrees of freedom per node.

    Returns
    -------
    ndarray or sparse matrix
        Updated global stiffness matrix.

    Algorithm
    ---------
    1. Compute the global DOF indices $I_e$ from the element topology $T_e$.
    2. Add the local stiffness matrix $K_e$ to the global stiffness matrix $K$ using index slicing $K[I_e, I_e] \mathrel{+}= K_e$.
    """
    element_nodes = topology_nodes(Te)
    indices = node_dof_indices(element_nodes, dof)
    element_matrix = as_float_array(Ke)
    if is_sparse(K):
        K[np.ix_(indices, indices)] = K[np.ix_(indices, indices)] + element_matrix
    else:
        K[np.ix_(indices, indices)] += element_matrix
    return K


def assmq(q, qe, Te, dof: int = 1):
    """
    Assemble one element force vector into the global internal-force vector.

    Parameters
    ----------
    q:
        Global internal-force vector updated in place.
    qe:
        Element internal-force vector.
    Te:
        Legacy FemLab topology row whose last entry stores the property id.
    dof:
        Degrees of freedom per node.

    Returns
    -------
    ndarray
        Updated global internal-force vector.

    Algorithm
    ---------
    1. Compute the global DOF indices $I_e$ from the element topology $T_e$.
    2. Reshape the element internal-force vector $q_e$ to a column vector.
    3. Add the element internal-force vector $q_e$ to the global internal-force vector $q$ using index slicing $q[I_e, 0] \mathrel{+}= q_e$.
    """
    element_nodes = topology_nodes(Te)
    indices = node_dof_indices(element_nodes, dof)
    q = as_float_array(q)
    qe = as_float_array(qe).reshape(-1, 1)
    q[indices, 0] += qe[:, 0]
    return q


__all__ = ["assmk", "assmq"]
