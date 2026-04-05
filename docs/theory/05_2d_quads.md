---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# 2D Isoparametric Quadrilateral Elements

This chapter explains the 4-node bilinear quadrilateral used throughout
`femlabpy`. The math is standard, but the focus here is on how the source code
organizes the parent-space derivatives, Jacobian solve, strain recovery, and
global assembly.

The implementation in `src/femlabpy/elements/quads.py`
uses a consistent set of array conventions:

- `Xe` is the element coordinate array with shape `(4, 2)` for a Q4 element.
- `Ge` stores the material row, usually `[E, nu]` or `[E, nu, type, t]`.
- Solid-mechanics displacement vectors are ordered `[u1, v1, u2, v2, u3, v3, u4, v4]`.
- Scalar potential kernels reuse the same geometry with one DOF per node.

## Bilinear shape functions

The parent domain is the square `[-1, 1] x [-1, 1]`. The four bilinear shape
functions are the tensor product of 1D linear Lagrange polynomials:

$$
N_1 = \frac{1}{4}(1-\xi)(1-\eta), \quad
N_2 = \frac{1}{4}(1+\xi)(1-\eta), \quad
N_3 = \frac{1}{4}(1+\xi)(1+\eta), \quad
N_4 = \frac{1}{4}(1-\xi)(1+\eta)
$$

These functions are used for both geometry and field interpolation. That is the
reason the element is called isoparametric.

The source computes derivatives in parent coordinates with `_q4_dN(r_i, r_j,
nnodes)`. For the public Q4 routines, `nnodes` is `4`, but the helper also
supports the derivative layout used by the serendipity Q8 branch in the same
module.

## Isoparametric mapping and the Jacobian

For a given Gauss point `(xi, eta)`, the physical coordinates and field values
are interpolated with the same shape functions:

$$
x(\xi, \eta) = \sum_i N_i x_i, \quad y(\xi, \eta) = \sum_i N_i y_i
$$

The source evaluates the Jacobian transpose directly as

```python
dN = _q4_dN(r[i], r[j], nnodes)
Jt = dN @ Xe
```

Here `dN` has shape `(2, 4)`:

- row 0 contains `dN/dxi`
- row 1 contains `dN/deta`

and `Jt` is a `2 x 2` matrix that maps parent-space derivatives to physical
derivatives.

The implementation intentionally uses `np.linalg.solve(Jt, dN)` instead of
`np.linalg.inv(Jt) @ dN`. That is the right choice for two reasons:

- it avoids forming an explicit inverse
- it is numerically safer when the element is distorted

The determinant of `Jt` is the area scaling factor that converts parent-space
quadrature weights into physical area weights.

## Strain-displacement matrix

For plane stress and plane strain, the strain vector is

$$
\boldsymbol{\varepsilon} = [\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]^T
$$

and the element displacement vector is

$$
\mathbf{d}_e = [u_1, v_1, u_2, v_2, u_3, v_3, u_4, v_4]^T
$$

The source builds the `B` matrix with a small loop over nodes:

```python
def _q4_B(dN):
    nnodes = dN.shape[1]
    B = np.zeros((3, 2 * nnodes), dtype=float)
    for k in range(nnodes):
        B[0, 2 * k] = dN[0, k]
        B[1, 2 * k + 1] = dN[1, k]
        B[2, 2 * k] = dN[1, k]
        B[2, 2 * k + 1] = dN[0, k]
    return B
```

That layout is deliberate:

- row 0 maps `du/dx`
- row 1 maps `dv/dy`
- row 2 maps the engineering shear strain `du/dy + dv/dx`

The same helper is reused in the scalar and elastoplastic kernels, so the
element math stays in one place.

## Gauss-Legendre integration

The Q4 element uses the standard `2 x 2` Gauss rule. The source stores the
abscissas and weights in `_q4_gauss_points()`:

```python
r = np.array([-1.0, 1.0], dtype=float) / np.sqrt(3.0)
w = np.array([1.0, 1.0], dtype=float)
```

Each of the four integration points is visited in a nested loop. The helper
`_q4_gp_index(i, j)` flattens the tensor-product indexing into the row order
used by the returned Gauss-point tables.

For a linear elastic problem, the stiffness matrix is accumulated as

$$
\mathbf{K}^e = \sum w_i w_j \mathbf{B}^T \mathbf{D} \mathbf{B} |\mathbf{J}| t
$$

where `t` is the thickness stored in `Ge[3]` when provided, otherwise `1.0`.

## Implementation in femlabpy

The public Q4 routines follow the same pattern:

- compute parent derivatives
- map them to physical derivatives with `np.linalg.solve`
- build `B`
- evaluate strains and stresses
- accumulate element contributions

### `keq4e`

The stiffness routine keeps the implementation close to the math:

```python
Xe = as_float_array(Xe)
props = as_float_array(Ge).reshape(-1)
plane_strain = props.size > 2 and int(props[2]) == 2
D = _plane_elastic_matrix_2d(props, plane_strain=plane_strain)
t = props[3] if props.size > 3 else 1.0
nnodes = Xe.shape[0]
r, w = _q4_gauss_points()
Ke = np.zeros((2 * nnodes, 2 * nnodes), dtype=float)
for i in range(2):
    for j in range(2):
        dN = _q4_dN(r[i], r[j], nnodes)
        Jt = dN @ Xe
        dN_global = np.linalg.solve(Jt, dN)
        B = _q4_B(dN_global)
        Ke += w[i] * w[j] * t * (B.T @ D @ B) * np.linalg.det(Jt)
```

The useful point is that the same code path works for plane
stress and plane strain. Only the constitutive matrix `D` changes.

### `qeq4e`

The stress recovery routine follows the same Gauss loop, but it stores the
element history in tables:

- `Se` has shape `(4, 3)` and stores `[sxx, syy, txy]`
- `Ee` has shape `(4, 3)` and stores `[exx, eyy, gxy]`
- `qe` has shape `(8, 1)` and stores the equivalent nodal internal force

The code updates the Gauss-point arrays one point at a time:

```python
gp = _q4_gp_index(i, j)
dN = _q4_dN(r[i], r[j], nnodes)
Jt = dN @ Xe
dN_global = np.linalg.solve(Jt, dN)
B = _q4_B(dN_global)
Ee[gp] = (B @ Ue).ravel()
Se[gp] = Ee[gp] @ D
qe += w[i] * w[j] * t * (B.T @ Se[gp].reshape(-1, 1)) * np.linalg.det(Jt)
```

That structure is useful because the Gauss-point values are kept in element
order, so the caller can inspect them without reconstructing any connectivity.

### Global assembly

The global wrappers `kq4e` and `qq4e` are thin loops around the element
kernels. Each row of `T` is treated as `[n1, n2, n3, n4, mat_id]`, and the
material row is selected with `topology_property(row)` before the local matrix
is assembled with `assmk` or `assmq`.

This makes the data flow explicit:

- `kq4e` builds the global stiffness matrix
- `qq4e` builds the global internal force vector and stores element stresses
- the same topology and material conventions are used in both

## Scalar potential version

The scalar kernels `keq4p`, `qeq4p`, `kq4p`, and `qq4p` use the same geometry
and the same Gauss integration strategy, but with one DOF per node.

The only real differences are:

- the gradient operator is `2 x 4` instead of `3 x 8`
- the element conductivity matrix is `4 x 4`
- the optional reaction coefficient `b` adds a consistent reaction term in
  `keq4p`

This is why the scalar and solid implementations can live in the same module:
they share the geometry and quadrature, but differ only in the constitutive and
DOF layout.

## Summary

The Q4 element in `femlabpy` is a good example of the library's implementation
style:

- keep the low-level geometry in small helpers
- evaluate derivatives with `np.linalg.solve`
- reuse the same Gauss loop for stiffness and recovery
- keep global assembly separate from element math

That split keeps the code readable and makes it easier to reason about the
exact shape of each array at every stage.
