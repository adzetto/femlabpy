# Chapter 1: Fundamentals & Data Structures

## 1.1 The femlabpy Philosophy

`femlabpy` is a pure-Python finite element library built around flat `numpy`
arrays rather than a deep object hierarchy. That choice is intentional. The
library follows the FemLab MATLAB/Scilab style where the model is represented by
tables and vectors, and each solver step maps directly onto array operations.

That design has three practical consequences:

1. The code mirrors the math. A stiffness matrix is a 2D array, a displacement
   field is a column vector, and a mesh is a small set of integer tables.
2. The data flow is easy to inspect. You can print the arrays before and after
   assembly or boundary-condition application and see exactly what changed.
3. The implementation stays close to NumPy and SciPy primitives, which keeps the
   kernels compact and makes vectorized assembly straightforward.

The public API exposed through `src/femlabpy/__init__.py`
collects the core pieces of that workflow at package level. Functions such as
`init`, `setload`, `setbc`, `assmk`, `assmq`, and the element kernels are all
imported there so users can work in the same style as the legacy FemLab scripts.

## 1.2 Core Data Structures

A finite element model in `femlabpy` is defined by a small set of arrays. The
library does not hide these arrays behind wrappers; instead, the arrays are the
model.

### 1.2.1 Node Coordinates Matrix

The matrix `X` stores the coordinates of every node in the global coordinate
system. For a 2D model with `N` nodes, `X` has shape `(N, 2)`. For a 3D model,
the shape is `(N, 3)`.

Each row is one node:

$$
\mathbf{X} =
\begin{bmatrix}
x_1 & y_1 \\
x_2 & y_2 \\
\vdots & \vdots \\
x_N & y_N
\end{bmatrix}
$$

In the code, `X[i]` means the coordinates of node `i + 1` in the legacy
one-based numbering used by the mesh tables. The Python array is still
zero-based, so any user-facing node number must be shifted down before it is
used as an index.

### 1.2.2 Topology Matrix

The topology matrix `T` defines element connectivity. Each row corresponds to
one element. The first columns list the one-based node numbers that belong to
the element, and the final column stores the property number used to look up the
material row in `G`.

For a Q4 mesh, `T` typically has shape `(E, 5)`:

$$
\mathbf{T} =
\begin{bmatrix}
n_{1,1} & n_{1,2} & n_{1,3} & n_{1,4} & p_1 \\
\vdots & \vdots & \vdots & \vdots & \vdots \\
n_{E,1} & n_{E,2} & n_{E,3} & n_{E,4} & p_E
\end{bmatrix}
$$

The topology row is the entry point for almost every element kernel. In the
helpers layer, `topology_nodes()` returns the node portion of the row and
`topology_property()` returns the material/property id. That split matters
because the last column is not part of the connectivity itself.

### 1.2.3 Material Properties

The material matrix `G` stores the physical properties associated with each
property id. Each row is one material definition. For common 2D elastic
problems, a row often looks like:

`[E, nu, type, t, rho]`

where:

`E`
: Young's modulus.

`nu`
: Poisson ratio.

`type`
: Material mode selector used by the element kernel, for example plane stress
  or plane strain.

`t`
: Thickness for 2D formulations.

`rho`
: Density used by dynamic mass-matrix kernels.

The helper `material_row()` converts the one-based property number stored in
`T` into the correct zero-based row of `G`.

### 1.2.4 Boundary Constraints

The boundary-condition table `C` stores prescribed displacements. Each row
identifies a node, a local degree of freedom, and the prescribed value.

For a 2D problem, the standard row format is:

`[node_id, dof_id, prescribed_value]`

The DOF numbering is one-based in the public tables:

`1`
: First translational DOF, usually `u_x`.

`2`
: Second translational DOF, usually `u_y`.

`3`
: Third translational DOF for 3D problems, usually `u_z`.

The boundary solver uses these rows to compute the flattened global DOF index
before it modifies the system.

### 1.2.5 Point Loads

The point-load table `P` has the same indexing convention as `C`, but the third
column stores a force component rather than a prescribed displacement.

`[node_id, dof_id, force_value]`

This table is consumed by `setload()`, which maps each row into the correct
entry of the global load vector.

### 1.2.6 Global Vectors and Matrices

The array-based workflow only becomes useful once the global system is defined.
`core.init()` allocates the containers used by the solvers:

- `K`: global stiffness matrix.
- `M`: global mass matrix, when dynamic analysis is requested.
- `p`: global external load vector.
- `q`: global internal force vector.

The helper returns dense `numpy.ndarray` objects by default, but it can switch
to sparse storage for large problems. The purpose is simple: the rest of the
library should not care whether the system is dense or sparse, only that the
container has the right shape and supports incremental assembly.

## 1.3 Array Helpers and Index Flow

The small helper functions in `src/femlabpy/_helpers.py`
do most of the indexing work. They are the bridge between the human-facing
tables and the zero-based NumPy arrays used by the solvers.

`as_float_array()` and `as_int_array()` normalize incoming data so the kernels do
not need to repeat type checks. `as_column()` forces a column-vector shape, which
is important because many solver routines expect `(ndof, 1)` rather than a flat
vector.

`rows()` and `cols()` provide MATLAB-like shape queries. They are useful when the
library accepts both row vectors and column vectors but still needs to preserve
the classroom conventions of the legacy FemLab scripts.

The most important indexing helpers are these:

- `element_dof_indices(node_numbers, dof)` expands one-based node numbers into
  flattened zero-based DOF indices.
- `node_dof_indices(node_numbers, dof)` is the flattened version used by
  assembly and load application.
- `topology_nodes(topology_row)` strips the property id from a legacy topology
  row.
- `material_row(materials, property_number)` resolves a one-based property id to
  the correct material row.

That flow is the same everywhere in the package:

1. Read a one-based user table.
2. Convert it to zero-based integer indices.
3. Use NumPy advanced indexing or sparse assembly to update the global arrays.

## 1.4 The Finite Element Analysis Sequence

The standard static workflow in `femlabpy` is intentionally explicit:

1. Allocate global arrays with `init()`.
2. Compute element matrices and assemble them into the global system.
3. Map the load table into the global right-hand side.
4. Apply constraints.
5. Solve the linear system.
6. Recover internal forces, stresses, and strains.

That sequence is easier to understand if you look at the actual helper calls:

```python
import numpy as np
import femlabpy as fp

K, p, q = fp.init(nn, dof)
K = fp.kq4e(K, T, X, G)
p = fp.setload(p, P)
K_bc, p_bc, _ = fp.setbc(K, p, C, dof)
u = np.linalg.solve(K_bc, p_bc)
q, S, E = fp.qq4e(np.zeros_like(p), T, X, G, u)
```

The important point is not the syntax. It is the data flow:

- `init()` creates a system of the correct size.
- `kq4e()` reads `X`, `T`, and `G`, builds one element matrix at a time, and
  scatters those contributions into `K`.
- `setload()` writes the forces into the right rows of `p`.
- `setbc()` modifies `K` and `p` in place so the prescribed values are enforced.
- `qq4e()` evaluates the element response after the displacement vector is known.

## 1.5 Worked Example

This example shows the same data flow using a single Q4 element. The goal is not
to solve a large model. The goal is to show exactly how the arrays are shaped
and how the indexing works.

### 1.5.1 Build the tables

```{code-block} python
import numpy as np
import femlabpy as fp

X = np.array([
    [0.0, 0.0],
    [2.0, 0.0],
    [2.0, 2.0],
    [0.0, 2.0],
])

T = np.array([
    [1, 2, 3, 4, 1],
])

G = np.array([
    [200e9, 0.3, 1, 0.1, 7850.0],
])

C = np.array([
    [1, 1, 0.0],
    [1, 2, 0.0],
    [4, 1, 0.0],
    [4, 2, 0.0],
])

P = np.array([
    [3, 2, -10000.0],
])
```

`X` stores coordinates row by row. `T` says that the single element uses nodes
1, 2, 3, and 4, and that it should read material row 1 from `G`. `C` and `P`
use the same one-based node numbering, which keeps the input tables readable and
compatible with the legacy FemLab conventions.

### 1.5.2 Read the element coordinates

```{code-block} python
element_nodes = T[0, :4]
element_coords = X[element_nodes - 1]
```

`T[0, :4]` extracts the four node ids of the first element. The subtraction by
`1` is the only translation between the user-facing table and the Python array.
Once the indices are zero-based, `X[element_nodes - 1]` returns the coordinates
needed by the element kernel.

### 1.5.3 Assemble the global system

```{code-block} python
nn = X.shape[0]
dof = 2
K, p, q = fp.init(nn, dof)

K = fp.kq4e(K, T, X, G)
p = fp.setload(p, P)
K_bc, p_bc, ks = fp.setbc(K, p, C, dof)
u = np.linalg.solve(K_bc, p_bc)
q, S, E = fp.qq4e(np.zeros_like(p), T, X, G, u)
```

`init()` creates an 8-by-8 stiffness matrix and 8-by-1 vectors because there are
4 nodes and 2 DOFs per node. `setload()` writes the downward force into the
global load vector using the flattened index

$$
(node - 1) \times dof + (local\_dof - 1)
$$

so the force at node 3 in the `y` direction lands in the correct slot. `setbc()`
then modifies the global system so the constrained DOFs behave as fixed
supports. Finally, `qq4e()` recovers the internal response from the solved
displacement field.

### 1.5.4 Inspect the flattened indices

```{code-block} python
from femlabpy._helpers import node_dof_indices

global_dofs = node_dof_indices([1, 4], dof=2)
```

For nodes 1 and 4 in a 2D model, `global_dofs` becomes `[0, 1, 6, 7]`. That is
the exact index set used by the low-level assembly and constraint routines. The
important pattern is that node numbering stays one-based in the input tables,
but every time the library touches a NumPy array it converts the table to
zero-based flattened DOFs.

## 1.6 What To Read Next

Once the array model is clear, the next chapters build on it directly:

- Chapter 2 explains the element kernels and how their local matrices are
  formed.
- Chapter 3 shows how local matrices become global systems through assembly.
- Chapter 5 extends the same data structures to dynamic analysis.
