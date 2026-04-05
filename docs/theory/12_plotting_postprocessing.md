# Plotting and Postprocessing

This chapter explains how `femlabpy` turns solved arrays into figures and how it
extracts support reactions from the internal force vector. The code is small, but
the conventions matter.

## Mesh and Field Plotting

The plotting helpers in `plotting.py` work directly on the topology and
coordinate arrays used elsewhere in the package.

### `plotelem`

`plotelem(T, X, ...)` draws the undeformed mesh. It reads each element row,
subtracts one from the one-based node ids, and plots the resulting polygon or
polyline.

The helper also supports optional node and element labels. That is useful when
you are checking topology ordering or debugging a bad mesh import.

### `plotforces`

`plotforces(T, X, P, ...)` draws nodal loads as arrows. The element table is not
used directly; it is accepted only for a consistent plotting signature.

The arrow size is scaled from the mesh span so a small model and a large model
remain readable with the same function call.

### `plotbc`

`plotbc(T, X, C, ...)` visualizes prescribed displacements. Zero-valued
constraints are shown as markers, while nonzero prescribed values are drawn as
arrows.

That distinction matches the actual solver behavior:

- zero values are fixed degrees of freedom;
- nonzero values are enforced displacements.

### `plotu`

`plotu(T, X, u, ...)` draws a scalar nodal field over a 2D or 3D mesh. In 2D it
uses `PolyCollection`. In 3D it uses `Poly3DCollection`.

The function does not attempt to infer the meaning of the scalar field. It just
colors each element by the mean value of its nodal samples.

## Element Contours

### `plotq4`

`plotq4(T, X, S, scomp, ...)` reconstructs a nodal contour from Q4 Gauss-point
results. The code expects the Q4 storage layout used by the Q4 recovery
functions: each component is stored at four Gauss points per element.

The implementation does three things:

1. it extracts the requested component with one-based `scomp`;
2. it maps the four Gauss-point values to corner values with a small
   extrapolation matrix;
3. it averages contributions from adjacent elements at each node.

The final contour is displayed with `tripcolor(..., shading="gouraud")` so the
result looks smooth across the mesh.

### `plott3`

`plott3(T, X, S, scomp, ...)` follows the same idea for T3 output. The recovery
is simpler because the T3 storage is already element-centered in the form used
by the rest of the codebase.

## Reactions

The `reaction()` function in `postprocess.py` reads the internal force vector
and extracts the support reactions associated with constrained degrees of
freedom.

It does not recompute the finite element solution. It assumes you already have
the assembled internal force vector `q` and the constraint table `C`.

### Inputs and outputs

- `q` is the global internal force vector.
- `C` is the constraint table.
- `dof` tells the helper how to flatten node and local DOF ids.
- `comp` optionally filters the returned rows by a single constrained
  component.

The output is either:

- `[node, local_dof, reaction]`, or
- `[constraint_row, reaction]` when `comp` is supplied.

### Why this helper is small

Support reactions are a recovery task, not a solve task. The helper exists so
the solve path can stay focused on matrix assembly while the postprocess layer
handles the bookkeeping required to report reactions cleanly.

## Typical Workflow

```python
import numpy as np
import femlabpy as fp

K, p, q = fp.init(nn, dof)
K = fp.kq4e(K, T, X, G)
p = fp.setload(p, P)
K, p, _ = fp.setbc(K, p, C, dof)
u = np.linalg.solve(K, p)
q, S, E = fp.qq4e(q, T, X, G, u)
R = fp.reaction(q, C, dof)
```

Then use the plotting helpers to inspect the result:

- `plotelem` for geometry,
- `plotforces` for loads,
- `plotbc` for constraints,
- `plotu` for scalar nodal fields,
- `plotq4` or `plott3` for stress or strain contours.
