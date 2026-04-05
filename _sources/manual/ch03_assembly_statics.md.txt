# Chapter 3: Assembly and Constraints

This chapter explains how `femlabpy` turns element-level matrices and vectors into a global system, then applies loads and constraints before solving. The key idea is simple: each element works in its own local numbering, but the solver always operates on one global set of degrees of freedom.

## 3.1 Global Assembly

Finite element assembly is a scatter operation. Each element produces a local stiffness matrix `Ke` and, when needed, a local internal-force vector `qe`. Those local quantities must be mapped into the correct rows and columns of the global arrays.

In the classical matrix form, the assembled stiffness matrix can be written as:

$$ K = \sum_e L_e^T K_e L_e $$

In the actual Python implementation, `femlabpy` does not build explicit Boolean matrices for every element. Instead, it computes global degree-of-freedom indices directly and uses NumPy indexing to place the local blocks.

### 3.1.1 How the indices are formed

The numbering convention is:

- Node numbers in the legacy FemLab tables are one-based.
- Python arrays are zero-based.
- Each node owns `dof` consecutive global degrees of freedom.

The helper `node_dof_indices` in `src/femlabpy/_helpers.py` expands a list of node numbers into flat global DOF indices. For a 2D problem with `dof = 2`, node `n` maps to indices:

```python
[2*(n-1), 2*(n-1) + 1]
```

So a 4-node quadrilateral in 2D becomes an 8-entry index list. Those eight indices are then used to extract or update the corresponding `8 x 8` submatrix.

This is the pattern used throughout the library:

```python
indices = node_dof_indices(topology_nodes(Te), dof)
K[np.ix_(indices, indices)] += Ke
q[indices, 0] += qe[:, 0]
```

The `np.ix_` call is important. It converts a flat list of DOF indices into a two-dimensional open mesh so NumPy selects the full row/column block instead of pairing indices one by one.

### 3.1.2 What `assmk` and `assmq` do

The low-level assembly helpers are intentionally small:

- `assmk(K, Ke, Te, dof)` inserts one element stiffness matrix into the global stiffness matrix.
- `assmq(q, qe, Te, dof)` adds one element internal-force vector into the global internal-force vector.

The implementation in `src/femlabpy/assembly.py` does exactly two things:

1. Extract the node numbers from the topology row `Te`.
2. Convert those node numbers into global DOF indices and accumulate the local values at those positions.

That design keeps the element kernels simple. The element routine computes the local matrix or vector, and the assembly helper handles only the global scatter step.

### 3.1.3 Dense and sparse assembly

For small problems, a dense `numpy.ndarray` is acceptable. For realistic models, the global matrix is mostly zeros and should be stored sparsely.

`femlabpy` uses `scipy.sparse.lil_matrix` when a matrix is being built incrementally. LIL is a good construction format because repeated single-block insertions are cheap. After assembly, the matrix should be converted to `csr` or `csc` before solving.

```python
import numpy as np
import scipy.sparse as sp

def create_global_matrix(ndof, use_sparse=True):
    if use_sparse:
        return sp.lil_matrix((ndof, ndof), dtype=float)
    return np.zeros((ndof, ndof), dtype=float)

def add_element_matrix(K, Ke, edof):
    ix = np.ix_(edof, edof)
    K[ix] += Ke
    return K
```

The subtle point is that sparse and dense code can share the same indexing logic. Only the container type changes.

### 3.1.4 A small assembly example

This example mirrors the indexing pattern used in the library:

```python
import numpy as np

# Two elements, one DOF per node
K = np.zeros((3, 3))

ke1 = np.array([[1.0, -1.0],
                [-1.0, 1.0]])
ke2 = np.array([[2.0, -2.0],
                [-2.0, 2.0]])

edof1 = np.array([0, 1])
edof2 = np.array([1, 2])

K[np.ix_(edof1, edof1)] += ke1
K[np.ix_(edof2, edof2)] += ke2
```

The second element shares node 1 with the first element, so its contribution is added to the same global row and column. That overlap is the whole point of assembly.

## 3.2 Dirichlet Boundary Conditions

Dirichlet boundary conditions prescribe displacements. In structural problems, they are used to fix supports, remove rigid-body modes, or enforce known motions.

The implementation in `src/femlabpy/boundary.py` uses direct elimination with a replacement diagonal value `ks`. This is not a pure penalty method in the textbook sense. It also transfers the coupling terms from the constrained column into the right-hand side before zeroing the row and column, so nonzero prescribed displacements are handled correctly.

### 3.2.1 How `setbc` works

The function `setbc(K, p, C, dof)` expects a legacy constraint table:

- For `dof = 1`, each row is `[node, value]`.
- For `dof > 1`, each row is `[node, local_dof, value]`.

The code does the following for each constrained DOF:

1. Compute the global DOF index from the node number and local DOF.
2. Compute `ks = 0.1 * max_abs_diagonal(K)`.
3. Subtract the column contribution from the load vector when the prescribed value is nonzero.
4. Zero the corresponding row and column.
5. Put `ks` on the diagonal and set the load entry to `ks * value`.

That means the constrained equation becomes `ks * u_i = ks * value`, so the solution is forced to the prescribed displacement while the rest of the system sees the corrected coupling forces.

```python
def enforce_dirichlet(K, p, dof_index, value, ks):
    if value != 0.0:
        p[:, 0] -= K[:, dof_index] * value
    K[dof_index, :] = 0.0
    K[:, dof_index] = 0.0
    K[dof_index, dof_index] = ks
    p[dof_index, 0] = ks * value
    return K, p
```

### 3.2.2 Why this matters

Without boundary enforcement, the global stiffness matrix is often singular because rigid-body modes are still free. Applying `setbc` removes those unconstrained motions from the numerical problem while preserving the effect of a nonzero support displacement.

The same mechanism is used whether the problem is scalar or vector-valued. Only the mapping from node number to global DOF changes.

## 3.3 General Linear Constraints

Some constraints are not simple fixed displacements. Periodic boundary conditions, rigid links, and tying relations are all linear constraints of the form:

$$ G u = Q $$

Here `G` is the constraint matrix and `Q` is the constraint right-hand side. The solver in `src/femlabpy/boundary.py` forms the augmented saddle-point system:

$$
\begin{bmatrix}
K & G^T \\
G & 0
\end{bmatrix}
\begin{Bmatrix}
u \\
\lambda
\end{Bmatrix}
=
\begin{Bmatrix}
p \\
Q
\end{Bmatrix}
$$

The unknown `lambda` contains the constraint forces.

### 3.3.1 What `solve_lag_general` adds

`solve_lag_general(K, p, G, Q)` is the general solver behind the higher-level constrained routines. It does three practical things before solving:

1. It scales the constraint rows by a factor derived from the largest diagonal entry of `K`.
2. It assembles the augmented matrix with either `np.block` for dense systems or `scipy.sparse.bmat` for sparse systems.
3. It solves the enlarged linear system and optionally returns the physical Lagrange multipliers.

That scaling step is easy to miss, but it matters. Constraint rows can be numerically tiny compared with stiffness rows, and the scaling keeps the block system better conditioned.

### 3.3.2 Practical interpretation

If `G u = Q` encodes a rigid tie, the solver enforces the tie exactly. If `G` encodes periodic boundary relations, the same machinery can be used to relate opposing boundaries without manually eliminating DOFs.

In other words, the direct elimination route is for fixed supports, and the augmented route is for linear coupling between unknowns.

## 3.4 Loads and Reactions

Loads are stored in a global vector `p`, but the input tables in `femlabpy` are node-based. The helper functions in `src/femlabpy/loads.py` convert those node tables into DOF-level updates.

### 3.4.1 `setload`

`setload(p, P)` replaces the load values at the listed nodes. This is the right choice when the table fully defines the applied load pattern.

The code computes the flattened DOF indices from the node numbers and writes the listed values directly into the global vector.

```python
indices = ((loads[:, [0]].astype(int) - 1) * dof + np.arange(dof)).reshape(-1)
p[indices, 0] = loads[:, 1:1 + dof].reshape(-1)
```

This is a one-to-one mapping: whatever appears in `P` becomes the new value at those DOFs.

### 3.4.2 `addload`

`addload(p, P)` accumulates loads instead of replacing them. That matters when several independent load cases contribute to the same node or when loads are assembled in stages.

The implementation uses `np.add.at`, which is the safe choice when the same global DOF appears more than once in the list.

### 3.4.3 Reactions

After solving, support reactions are extracted from the global internal-force vector using `reaction(q, C, dof)` in `src/femlabpy/postprocess.py`.

The function maps the constrained entries back to the original boundary-condition table and returns a compact reaction list. In the vector-valued case, the output columns are `[node, local_dof, reaction]`.

That means the solution process is:

1. Assemble `K` and `p`.
2. Apply boundary conditions.
3. Solve for `u`.
4. Recover reactions from the internal-force vector at constrained DOFs.

## 3.5 Full example

This minimal example shows the same mechanics used by the library: element scatter, support enforcement, and reaction recovery.

```python
import numpy as np

# One-dimensional, three-node system
K = np.zeros((3, 3), dtype=float)
p = np.zeros((3, 1), dtype=float)

# Element matrices
ke1 = np.array([[1.0, -1.0], [-1.0, 1.0]])
ke2 = np.array([[2.0, -2.0], [-2.0, 2.0]])

# Assemble
K[np.ix_([0, 1], [0, 1])] += ke1
K[np.ix_([1, 2], [1, 2])] += ke2

# Apply a nodal load at DOF 2
p[1, 0] = 100.0

# Keep the assembled system before boundary conditions for reaction recovery
K_before_bc = K.copy()
p_before_bc = p.copy()

# Constrain DOF 0 to zero
ks = 0.1 * np.max(np.abs(np.diag(K)))
K[0, :] = 0.0
K[:, 0] = 0.0
K[0, 0] = ks
p[0, 0] = 0.0

u = np.linalg.solve(K, p)
reaction_at_support = (K_before_bc @ u - p_before_bc)[0, 0]
```

The important part is not the toy numbers. It is the sequence: assemble, apply loads, constrain, solve, then recover the reactions from the solved state.
