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

```{figure} ../_static/figures/fig_t3_element.png
:align: center
:width: 480px

**Figure 4.1:** 3-node constant strain triangle (T3/CST) with counter-clockwise node ordering and 2 DOFs per node.
```

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

### Area (barycentric) coordinates

The triangle uses area coordinates $L_1, L_2, L_3$. For a point $P$ inside a triangle with vertices $1, 2, 3$:

$$
L_i = \frac{A_i}{A}
$$

where $A$ is the total area and $A_i$ is the area of the sub-triangle formed by replacing vertex $i$ with point $P$.

By construction, these coordinates satisfy:

$$
L_1 + L_2 + L_3 = 1
$$

### Derivation of shape functions

The position of any point in the triangle is interpolated as:

$$
\mathbf{x} = L_1 \mathbf{x}_1 + L_2 \mathbf{x}_2 + L_3 \mathbf{x}_3
$$

Expanding in $x$ and $y$ components with the partition of unity:

$$
\begin{cases}
x = L_1 x_1 + L_2 x_2 + (1 - L_1 - L_2) x_3 \\
y = L_1 y_1 + L_2 y_2 + (1 - L_1 - L_2) y_3
\end{cases}
$$

Rearranging:

$$
\begin{cases}
x - x_3 = L_1 (x_1 - x_3) + L_2 (x_2 - x_3) \\
y - y_3 = L_1 (y_1 - y_3) + L_2 (y_2 - y_3)
\end{cases}
$$

In matrix form:

$$
\begin{bmatrix} x - x_3 \\ y - y_3 \end{bmatrix} = \begin{bmatrix} x_1 - x_3 & x_2 - x_3 \\ y_1 - y_3 & y_2 - y_3 \end{bmatrix} \begin{bmatrix} L_1 \\ L_2 \end{bmatrix}
$$

Inverting to solve for $L_1, L_2$:

$$
\begin{bmatrix} L_1 \\ L_2 \end{bmatrix} = \frac{1}{2A} \begin{bmatrix} y_2 - y_3 & x_3 - x_2 \\ y_3 - y_1 & x_1 - x_3 \end{bmatrix} \begin{bmatrix} x - x_3 \\ y - y_3 \end{bmatrix}
$$

where $2A = (x_1 - x_3)(y_2 - y_3) - (x_2 - x_3)(y_1 - y_3)$ is twice the triangle area.

### Shape function gradients

Since $N_i = L_i$, the gradients follow directly from the inverse mapping:

$$
\frac{\partial N_1}{\partial x} = \frac{y_2 - y_3}{2A}, \quad \frac{\partial N_1}{\partial y} = \frac{x_3 - x_2}{2A}
$$

$$
\frac{\partial N_2}{\partial x} = \frac{y_3 - y_1}{2A}, \quad \frac{\partial N_2}{\partial y} = \frac{x_1 - x_3}{2A}
$$

$$
\frac{\partial N_3}{\partial x} = \frac{y_1 - y_2}{2A}, \quad \frac{\partial N_3}{\partial y} = \frac{x_2 - x_1}{2A}
$$

These gradients are **constant** over the element because they depend only on the vertex coordinates, not on position. This is why the CST produces a constant strain field.

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

### Derivation from kinematics

For 2D plane problems, the strain-displacement relations are:

$$
\varepsilon_{xx} = \frac{\partial u}{\partial x}, \quad \varepsilon_{yy} = \frac{\partial v}{\partial y}, \quad \gamma_{xy} = \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}
$$

Using the finite element interpolation $u = \sum_i N_i u_i$ and $v = \sum_i N_i v_i$:

$$
\varepsilon_{xx} = \sum_i \frac{\partial N_i}{\partial x} u_i, \quad \varepsilon_{yy} = \sum_i \frac{\partial N_i}{\partial y} v_i
$$

$$
\gamma_{xy} = \sum_i \left( \frac{\partial N_i}{\partial y} u_i + \frac{\partial N_i}{\partial x} v_i \right)
$$

For the element displacement vector:

$$
\mathbf{u}^e = [u_1, v_1, u_2, v_2, u_3, v_3]^T
$$

The strain vector in Voigt form is:

$$
\boldsymbol{\varepsilon} = [\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]^T = \mathbf{B} \mathbf{u}^e
$$

### Strain-displacement matrix structure

Using the notation $b_i = y_{j} - y_{k}$ and $c_i = x_{k} - x_{j}$ (cyclic permutation), where indices $(i,j,k)$ cycle through $(1,2,3)$:

$$
\mathbf{B} = \frac{1}{2A}
\begin{bmatrix}
b_1 & 0 & b_2 & 0 & b_3 & 0 \\
0 & c_1 & 0 & c_2 & 0 & c_3 \\
c_1 & b_1 & c_2 & b_2 & c_3 & b_3
\end{bmatrix}
$$

This can be written as $\mathbf{B} = [\mathbf{B}_1 | \mathbf{B}_2 | \mathbf{B}_3]$ where each $\mathbf{B}_i$ is:

$$
\mathbf{B}_i = \frac{1}{2A} \begin{bmatrix} b_i & 0 \\ 0 & c_i \\ c_i & b_i \end{bmatrix}
$$

### Why the CST is constant strain

Since $b_i$ and $c_i$ depend only on nodal coordinates (not on position within the element), $\mathbf{B}$ is constant. Therefore:

$$
\boldsymbol{\varepsilon} = \mathbf{B} \mathbf{u}^e = \text{constant over the element}
$$

This is the defining characteristic of the Constant Strain Triangle.

## Element stiffness derivation

### Principle of virtual work

The element stiffness matrix is derived from the strain energy:

$$
U = \frac{1}{2} \int_{\Omega^e} \boldsymbol{\varepsilon}^T \mathbf{D} \boldsymbol{\varepsilon} \, dA = \frac{1}{2} (\mathbf{u}^e)^T \left( \int_{\Omega^e} \mathbf{B}^T \mathbf{D} \mathbf{B} \, dA \right) \mathbf{u}^e
$$

Since $\mathbf{B}$ is constant over the element:

$$
\mathbf{K}^e = \int_{\Omega^e} \mathbf{B}^T \mathbf{D} \mathbf{B} \, dA = \mathbf{B}^T \mathbf{D} \mathbf{B} \cdot A
$$

This is exact—no numerical integration is needed.

### Implementation in `ket3e`

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

---

## Appendix: Constitutive Matrix Derivations

### A.1 3D Hooke's Law

The full 3D stress-strain relationship for isotropic linear elasticity:

$$
\boldsymbol{\sigma} = \mathbf{D}^{3D} \boldsymbol{\varepsilon}
$$

In matrix form:

$$
\begin{bmatrix} \sigma_{xx} \\ \sigma_{yy} \\ \sigma_{zz} \\ \tau_{xy} \\ \tau_{yz} \\ \tau_{xz} \end{bmatrix}
= \frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix}
1-\nu & \nu & \nu & 0 & 0 & 0 \\
\nu & 1-\nu & \nu & 0 & 0 & 0 \\
\nu & \nu & 1-\nu & 0 & 0 & 0 \\
0 & 0 & 0 & \frac{1-2\nu}{2} & 0 & 0 \\
0 & 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 \\
0 & 0 & 0 & 0 & 0 & \frac{1-2\nu}{2}
\end{bmatrix}
\begin{bmatrix} \varepsilon_{xx} \\ \varepsilon_{yy} \\ \varepsilon_{zz} \\ \gamma_{xy} \\ \gamma_{yz} \\ \gamma_{xz} \end{bmatrix}
$$

### A.2 Plane stress derivation

**Assumption:** $\sigma_{zz} = \tau_{xz} = \tau_{yz} = 0$ (thin plate in $xy$-plane)

From the third row of 3D Hooke's law:

$$
\sigma_{zz} = \frac{E}{(1+\nu)(1-2\nu)} \left[ \nu\varepsilon_{xx} + \nu\varepsilon_{yy} + (1-\nu)\varepsilon_{zz} \right] = 0
$$

Solving for $\varepsilon_{zz}$:

$$
\varepsilon_{zz} = -\frac{\nu}{1-\nu}(\varepsilon_{xx} + \varepsilon_{yy})
$$

Substituting back into the $\sigma_{xx}$ equation:

$$
\sigma_{xx} = \frac{E}{(1+\nu)(1-2\nu)} \left[ (1-\nu)\varepsilon_{xx} + \nu\varepsilon_{yy} + \nu\varepsilon_{zz} \right]
$$

$$
= \frac{E}{(1+\nu)(1-2\nu)} \left[ (1-\nu)\varepsilon_{xx} + \nu\varepsilon_{yy} - \frac{\nu^2}{1-\nu}(\varepsilon_{xx} + \varepsilon_{yy}) \right]
$$

Simplifying (after algebra):

$$
\sigma_{xx} = \frac{E}{1-\nu^2} \left[ \varepsilon_{xx} + \nu\varepsilon_{yy} \right]
$$

The resulting **plane stress** constitutive matrix:

$$
\mathbf{D}_{ps} = \frac{E}{1-\nu^2}
\begin{bmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \frac{1-\nu}{2}
\end{bmatrix}
$$

### A.3 Plane strain derivation

**Assumption:** $\varepsilon_{zz} = \gamma_{xz} = \gamma_{yz} = 0$ (thick body constrained in $z$-direction)

Direct substitution of $\varepsilon_{zz} = 0$ into 3D Hooke's law:

$$
\sigma_{xx} = \frac{E}{(1+\nu)(1-2\nu)} \left[ (1-\nu)\varepsilon_{xx} + \nu\varepsilon_{yy} \right]
$$

$$
\sigma_{yy} = \frac{E}{(1+\nu)(1-2\nu)} \left[ \nu\varepsilon_{xx} + (1-\nu)\varepsilon_{yy} \right]
$$

$$
\tau_{xy} = \frac{E}{2(1+\nu)} \gamma_{xy} = G \gamma_{xy}
$$

The **plane strain** constitutive matrix:

$$
\mathbf{D}_{p\varepsilon} = \frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix}
1-\nu & \nu & 0 \\
\nu & 1-\nu & 0 \\
0 & 0 & \frac{1-2\nu}{2}
\end{bmatrix}
$$

### A.4 Out-of-plane stress in plane strain

In plane strain, $\sigma_{zz} \neq 0$:

$$
\sigma_{zz} = \nu(\sigma_{xx} + \sigma_{yy}) = \frac{E\nu}{(1+\nu)(1-2\nu)}(\varepsilon_{xx} + \varepsilon_{yy})
$$

This stress must be considered in yield criteria (e.g., von Mises).

### A.5 Relationship between plane stress and plane strain

The two formulations are related through the substitution:

| Plane Stress | Plane Strain Equivalent |
|--------------|-------------------------|
| $E$ | $\frac{E}{1-\nu^2}$ |
| $\nu$ | $\frac{\nu}{1-\nu}$ |

**Proof:** Starting from plane strain matrix with substitutions:

$$
\frac{E'}{(1+\nu')(1-2\nu')} = \frac{E/(1-\nu^2)}{(1+\frac{\nu}{1-\nu})(1-\frac{2\nu}{1-\nu})}
$$

Simplifying the denominator:

$$
\left(1+\frac{\nu}{1-\nu}\right) = \frac{1}{1-\nu}, \quad \left(1-\frac{2\nu}{1-\nu}\right) = \frac{1-3\nu}{1-\nu}
$$

After simplification, we recover the plane stress matrix. ∎
