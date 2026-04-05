# Chapter 6: Periodic Boundaries and Gmsh I/O

Periodic boundary conditions are the part of the workflow that turns a small unit cell into a usable representation of a repeating material. In `femlabpy`, that means three separate jobs have to work together:

1. Pair the nodes on opposite faces of the mesh.
2. Build constraint equations that express the macro strain.
3. Read meshes from Gmsh in a way that preserves node order, element tags, and explicit topology tables.

This chapter follows that flow from geometry checks to homogenization and mesh import.

## 6.1 Periodic boundary conditions

Periodic conditions are used when a representative volume element, or RVE, should deform like one tile in an infinite repeating pattern. The idea is simple: nodes on opposing boundaries must move in a compatible way, with the difference between the two sides driven by the imposed macro strain.

In `femlabpy`, the public API is axis-based. That matters because periodicity is not inferred from left and right boundary lists supplied by the caller. Instead, the library inspects the coordinates directly and pairs nodes along a chosen axis.

### Pairing nodes

The main entry point is `find_periodic_pairs(X, axis, tol=1e-6)`.

- `X` is the node coordinate array.
- `axis` selects the periodic direction: `0` for x, `1` for y, `2` for z.
- `tol` is scaled by the domain size, which keeps the matching rule stable when the mesh is large or very small.

The function finds nodes on the minimum and maximum boundary of the chosen axis, compares their transverse coordinates, and returns 1-based node pairs. If the mesh is 2D, the returned array has two columns: `[left_node, right_node]`.

That is the practical point to remember: the code does not need manually curated node lists. It needs a mesh that is actually aligned on opposite faces.

If you want to check more than one periodic axis, use `find_all_periodic_pairs(X, periodic_axes, tol=1e-6)`. In 2D, the common pattern is to request axes `[0, 1]` and then stack the returned pair arrays before building the constraint matrix.

### Mesh checks before solving

Before assembling constraints, call `check_periodic_mesh(X, axis, tol=1e-6)`.

The returned report tells you whether the opposite boundaries are structurally compatible:

- `valid` indicates whether the mesh passed the periodic check.
- `n_left` and `n_right` report how many nodes were found on each side.
- `max_mismatch` gives the worst pairing mismatch seen during validation.
- `message` explains the failure or confirms that pairing succeeded.

Use this as a fast guardrail. A mismatch in node counts usually means the meshing strategy is not periodic enough for an RVE solve, even if the geometry itself looks symmetric.

### Constraint matrix and macro strain

`periodic_constraints(X, pairs, dof, eps_macro=None)` converts paired nodes into the linear constraint system used by the solver.

The function builds two objects:

- `G`, the constraint matrix.
- `Q`, the right-hand-side vector.

For each pair and each degree of freedom, one row is created in `G`. The right node receives a `+1` and the left node receives a `-1`. If a macro strain is supplied, `Q` carries the displacement jump implied by that strain.

For 2D problems, the Voigt order is `[exx, eyy, gxy]`. Inside the implementation, the engineering shear strain is converted to tensor form by splitting it equally across the off-diagonal terms. That keeps the constraint consistent with the strain tensor used in the rest of the library.

`apply_macro_strain(X, pairs, eps_macro, dof)` is just a convenience wrapper when you only want `Q`.

### Stabilizing the rigid mode

Purely periodic systems still have a rigid translation unless you pin one corner or add an equivalent reference constraint. `fix_corner(X, C_existing, dof)` handles that by finding the node closest to the minimum coordinate corner and adding zero-displacement rows for the active DOFs.

That step is not cosmetic. Without it, the augmented system can remain singular even if the periodic equations are correct.

### Solving the periodic system

`solve_periodic(K, p, X, pairs, dof, eps_macro=None, return_lagrange=False)` wraps the constraint build and the augmented solve.

The solver is built around the block system

$$
\begin{bmatrix}
\mathbf{K} & \mathbf{G}^T \\
\mathbf{G} & \mathbf{0}
\end{bmatrix}
\begin{bmatrix}
\mathbf{u} \\
\boldsymbol{\lambda}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{p} \\
\mathbf{Q}
\end{bmatrix}
$$

In practice, `femlabpy` uses this for both direct periodic solves and homogenization. The Lagrange multipliers are useful if you want boundary reactions or if you want to inspect how strongly the constraints are being enforced.

### Homogenization workflow

`homogenize(K, T, X, G_mat, pairs, dof, element_type="q4")` computes the effective stiffness matrix of the RVE.

The function applies unit macro strain cases one by one, solves the constrained system, averages the resulting stress, and assembles the effective constitutive matrix column by column. For 2D, that matrix is `3 x 3`; for 3D, it is `6 x 6`.

The implementation is deliberately straightforward:

- Build a unit macro strain vector.
- Solve the periodic system for that strain.
- Compute the volume-averaged stress.
- Insert the result into one column of `C_eff`.

That is the right mental model for reading the code. The homogenizer is not doing anything mysterious. It is repeating one constrained solve per canonical strain state.

### Practical 2D workflow

The following pattern matches the current API and is the safest way to work with a periodic quadrilateral RVE:

```python
import numpy as np
import femlabpy as fp

mesh = fp.load_gmsh2("rve.msh")
X = mesh.positions[:, :2]
T = mesh.quads.astype(int)

report_x = fp.check_periodic_mesh(X, axis=0, tol=1e-6)
report_y = fp.check_periodic_mesh(X, axis=1, tol=1e-6)
if not report_x["valid"] or not report_y["valid"]:
    raise ValueError(report_x["message"] if not report_x["valid"] else report_y["message"])

pairs_x = fp.find_periodic_pairs(X, axis=0, tol=1e-6)
pairs_y = fp.find_periodic_pairs(X, axis=1, tol=1e-6)
pairs = np.vstack([pairs_x, pairs_y])

nn = mesh.nbNod
K, p, q = fp.init(nn, dof=2)

G = np.array([[210e9, 0.30, 1, 1.0]])
K = fp.kq4e(K, T, X, G)

eps_macro = np.array([0.01, 0.0, 0.0])
u = fp.solve_periodic(K, p, X, pairs, dof=2, eps_macro=eps_macro)
```

The same node pairs are reused for both constraint generation and homogenization. Once the pairing is correct, the rest of the workflow is just matrix assembly and solving.

## 6.2 Gmsh import workflow

`femlabpy` uses `io.gmsh` to convert `.msh` files into a normalized `GmshMesh` object. That object carries both the Python-friendly fields used internally and the legacy names expected by older FemLab examples.

### Choosing a loader

Two loaders are available:

- `load_gmsh(filename)` follows the legacy semantics and materializes all explicit element arrays.
- `load_gmsh2(filename, which=None)` is more flexible and lets you control which explicit topology arrays are created.

Use `load_gmsh2` when you want to keep the mesh object smaller. If `which` is `None`, all explicit arrays are loaded. If `which` is `-1` or an empty iterable, the explicit arrays are skipped.

### What the mesh object contains

The returned `GmshMesh` includes:

- `positions`, the node coordinates.
- `element_infos`, a compact per-element summary.
- `element_tags` and `element_nodes`, the generalized tag and node tables.
- `triangles`, `quads`, `tets`, `hexa`, and related explicit topology arrays.
- `bounds_min` and `bounds_max`, which are useful for quick domain checks.

The object also keeps legacy aliases such as `POS`, `TRIANGLES`, `QUADS`, `nbTriangles`, and similar names. That makes it easier to reuse older scripts while still using the modern Python fields.

### Physical tags and topology columns

The loader preserves Gmsh physical groups. For explicit topology arrays, the last column stores the first tag for that element. In other words, a quadrilateral row looks like this:

```text
[n1, n2, n3, n4, prop_id]
```

That convention is what lets the assembly routines map each element to the correct material row without extra bookkeeping.

### Version handling

The loader accepts legacy ASCII meshes directly. For modern Gmsh 4.x meshes, the code can convert through the official Gmsh SDK when the optional mesh extra is installed. If that package is missing, the loader raises a clear error instead of silently misreading the file.

That behavior is important. Periodic boundary workflows depend on exact node ordering and consistent element connectivity, so a partial parse would be worse than a hard failure.

### Practical import example

```python
import femlabpy as fp

mesh = fp.load_gmsh2("rve.msh")

print(mesh.nbNod)
print(mesh.nbElm)
print(mesh.quads.shape)
print(mesh.bounds_min, mesh.bounds_max)
```

If you only need a subset of explicit element tables, request them up front:

```python
mesh = fp.load_gmsh2("rve.msh", which=[3, 4])
```

That is useful when the mesh contains many element types but the analysis only uses quads and tets.

### Gmsh generation tips for periodic models

The loader can only preserve periodic compatibility if the mesh itself was built compatibly. For RVEs, the safest pattern in Gmsh is:

1. Build the geometry with opposite edges derived from the same parameterization.
2. Use a structured or transfinite mesh where possible.
3. Recombine when you want quads instead of triangles.
4. Check both periodic axes with `check_periodic_mesh` before solving.

If the boundary discretization is inconsistent, no amount of post-processing in Python will repair the pairing.

### End-to-end example

```python
import numpy as np
import femlabpy as fp

mesh = fp.load_gmsh2("unit_cell.msh")
X = mesh.positions[:, :2]
T = mesh.quads.astype(int)

left_right = fp.find_periodic_pairs(X, axis=0, tol=1e-6)
bottom_top = fp.find_periodic_pairs(X, axis=1, tol=1e-6)
pairs = np.vstack([left_right, bottom_top])

K, p, q = fp.init(mesh.nbNod, dof=2)
G = np.array([[70e9, 0.25, 1, 1.0]])
K = fp.kq4e(K, T, X, G)

C_eff = fp.homogenize(K, T, X, G, pairs, dof=2, element_type="q4")
print(C_eff)
```

This is the cleanest mental model for the whole chapter: verify the mesh, pair the boundaries, build the constraints, solve once for each macro strain state, and let the averaging routines assemble the effective response.
