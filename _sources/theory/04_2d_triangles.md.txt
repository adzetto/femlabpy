---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
jupytext_version: 1.14.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# 2D Triangular Elements

This chapter focuses on the 3-node constant strain triangle used by `femlabpy`.
The goal is not just to derive the element, but to show how the source code
organizes the geometry, the element matrices, and the global assembly path.

The implementation in `src/femlabpy/elements/triangles.py`
uses the same conventions throughout:

- Topology rows are 1-based and follow `[n1, n2, n3, mat_id]`.
- Element coordinates `Xe` are always passed as a `(3, 2)` array.
- Solid-mechanics displacements are ordered `[u1, v1, u2, v2, u3, v3]`.
- Scalar potential problems use one DOF per node and reuse the same geometry.

## Geometry and array layout

The CST is the simplest 2D solid element, but the code is still careful about
geometry because the element is only valid when the triangle has a positive
area and the node order is counter-clockwise.

For a triangle with vertices `(x1, y1)`, `(x2, y2)`, and `(x3, y3)`, the area
is computed from the edge differences. The source uses the helper
`_triangle_geometry(Xe)` for one element and `_triangle_batch_geometry(Xe)` for
many elements at once. Both helpers return the same information:

- edge-difference vectors used to build the gradients of the shape functions
- the absolute area of the triangle

That split matters because the element routines use the same geometry in three
different places:

- `ket3e` for one element stiffness matrix
- `kt3e` for batched global stiffness assembly
- `qet3e` and `qt3e` for strain/stress or gradient recovery

## Shape functions and gradients

The triangle uses area coordinates, also called barycentric coordinates.
Inside the element, the three shape functions are exactly those coordinates:

$$
N_1 = L_1, \quad N_2 = L_2, \quad N_3 = L_3
$$

The practical benefit is that the gradients are constant over the element.
That is why the CST produces a constant strain field.

The source computes those gradients from the edge-difference helpers rather than
by symbolic manipulation at runtime. For a single triangle, the derivative array
is built as:

```python
a, area = _triangle_geometry(Xe)
dN = (1.0 / (2.0 * area)) * np.column_stack([-a[:, 1], a[:, 0]]).T
```

The resulting `dN` has shape `(2, 3)`:

- row 0 is `dN/dx` for the three nodes
- row 1 is `dN/dy` for the three nodes

This layout is reused in both the solid and scalar-potential kernels.

## Strain-displacement matrix

For plane stress and plane strain, each node has two DOFs, so the local
displacement vector is

$$
\mathbf{u}^e = [u_1, v_1, u_2, v_2, u_3, v_3]^T
$$

The strain vector in Voigt form is

$$
\boldsymbol{\varepsilon} = [\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]^T
$$

The source assembles the `B` matrix directly from the derivative table instead
of branching on the three strain components. The layout is:

$$
\mathbf{B} = \frac{1}{2A}
\begin{bmatrix}
b_1 & 0 & b_2 & 0 & b_3 & 0 \\
0 & c_1 & 0 & c_2 & 0 & c_3 \\
c_1 & b_1 & c_2 & b_2 & c_3 & b_3
\end{bmatrix}
$$

In code, that becomes a small dense matrix multiply with no numerical
integration loop because the CST is exact with one constant `B` matrix.

## Element stiffness in `ket3e`

The single-element stiffness routine is the cleanest place to see the element
logic:

```python
a, area = _triangle_geometry(Xe)
dN = (1.0 / (2.0 * area)) * np.column_stack([-a[:, 1], a[:, 0]]).T
B = np.array(
    [
        [dN[0, 0], 0.0, dN[0, 1], 0.0, dN[0, 2], 0.0],
        [0.0, dN[1, 0], 0.0, dN[1, 1], 0.0, dN[1, 2]],
        [dN[1, 0], dN[0, 0], dN[1, 1], dN[0, 1], dN[1, 2], dN[0, 2]],
    ],
    dtype=float,
)
D = _elastic_matrix(props, plane_strain=plane_strain)
Ke = (B.T @ D @ B) * area
```

The important point is that `ket3e` never approximates the integral with Gauss
points. The CST has constant strain, so the exact element matrix is simply

$$
\mathbf{K}^e = A \mathbf{B}^T \mathbf{D} \mathbf{B}
$$

The optional third entry in `Ge` switches between plane stress and plane strain.
The function reads `Ge[2]` only when it is present and equal to `2`.

## Batched assembly in `kt3e`

The global assembler is where the implementation becomes more performance
oriented. The source avoids a Python loop over node pairs and instead builds the
element data for the full topology batch first:

```python
nodes = topology[:, :3].astype(int) - 1
materials = as_float_array(G)[topology[:, -1].astype(int) - 1]
edges, area = _triangle_batch_geometry(coordinates[nodes])
```

From there, the code constructs a batched `B` tensor with shape `(nel, 3, 6)`
and a batched constitutive tensor with shape `(nel, 3, 3)`. The actual matrix
product is done with `np.einsum`:

```python
element_matrices = area[:, None, None] * np.einsum(
    "eik,ekl,elj->eij", B.transpose(0, 2, 1), D, B
)
```

That line is the key step in the implementation. It means:

- `e` indexes the element number
- `i` and `j` are the local DOF indices
- `k` and `l` are the contracted strain-space indices

So instead of looping over elements in Python, the code performs the full batch
of `B^T D B` products in one vectorized call.

The scatter step uses `element_dof_indices(nodes, 2, one_based=False)` and then
either:

- `np.add.at` for dense global matrices
- a COO sparse accumulation path for sparse global matrices

That separation keeps the assembly code fast without changing the public API.

## Stress and flux recovery

The recovery routines follow the same geometry logic, but they return more than
just the stiffness matrix:

- `qet3e` computes one element's strain, stress, and internal force vector
- `qt3e` loops over all triangles and scatters the internal force contributions

The returned arrays are shaped so downstream code can treat each element as a
small table:

- `Se` is `(nel, 3)` for `[sxx, syy, txy]`
- `Ee` is `(nel, 3)` for `[exx, eyy, gxy]`

The batched recovery path mirrors the stiffness assembly path closely:

```python
E = np.einsum("eij,ej->ei", B, element_displacements)
S = np.einsum("ei,eij->ej", E, D)
element_vectors = area[:, None] * np.einsum("eij,ej->ei", B.transpose(0, 2, 1), S)
```

The arrays are intentionally kept in element-major order so the caller can store
or postprocess them without reconstructing the mesh connectivity.

## Scalar potential version

The same geometry is reused for scalar potential and heat-flow problems through
`ket3p`, `qet3p`, `kt3p`, and `qt3p`.

The difference is only the DOF count and the size of the derivative operator:

- solid mechanics uses a `3 x 6` `B` matrix
- scalar problems use a `2 x 3` gradient operator

The scalar stiffness matrix follows the same pattern:

$$
\mathbf{K}^e = A \mathbf{B}^T \mathbf{D} \mathbf{B}
$$

with `D = k I` for isotropic conductivity. An optional reaction term can be
added in `ket3p`, which is useful for Poisson-type problems.

## Summary

The triangle kernels in `femlabpy` are small, but they show the same design
pattern used throughout the library:

- compute geometry once
- build compact element tensors
- use dense linear algebra on small blocks
- scatter the results into global arrays with explicit indexing

That structure makes the CST easy to audit and also makes the batched assembly
path fast enough for larger meshes.
