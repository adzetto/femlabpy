---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# 3D Solid Elements: Tetrahedrons and Hexahedrons

This chapter follows the actual implementation in `femlabpy.elements.solids`. The code path is the same for every solid element: build shape-function derivatives, map them through the Jacobian, assemble the strain-displacement matrix, evaluate the elastic matrix, and then either return the element stiffness or recover stress and internal force. The difference between T4 and H8 is mostly how the derivatives are obtained and whether the integration is exact or Gauss based.

## 1. The 4-Node Tetrahedron Element

The T4 element is the simplest 3D solid in the library. Each element has 4 nodes, 3 translational degrees of freedom per node, and therefore 12 local degrees of freedom. The implementation assumes the usual topology row layout:

`[n1, n2, n3, n4, mat_id]`

The coordinate block passed to the kernel is the 4 x 3 array `Xe`, and the material row `Ge` is expected to contain the isotropic elastic properties used by `_elastic3d_matrix`.

### 1.1 Constant derivatives and data flow

For T4, the parent-space derivatives are constant. That is why the code can define

```python
dN = np.array(
    [[1.0, 0.0, 0.0, -1.0],
     [0.0, 1.0, 0.0, -1.0],
     [0.0, 0.0, 1.0, -1.0]],
    dtype=float,
)
```

and reuse it directly in every call. The Jacobian is then

```python
J = dN @ Xe
```

which means the geometry of the element is carried entirely by `Xe`. Once `J` is known, the code pushes the derivatives into global coordinates with

```python
dN = np.linalg.solve(J, dN)
```

That line is the key state transition in the element: local derivatives become physical derivatives. From there, `_solid_B(dN)` builds the 6 x 12 strain-displacement matrix and `_elastic3d_matrix(Ge)` builds the 6 x 6 constitutive matrix.

### 1.2 Stiffness and internal force

`keT4e` computes the stiffness directly from the exact volume integral. In code, the integration reduces to

```python
return 2.0 * (B.T @ D @ B) * np.linalg.det(J)
```

The element is constant strain, so there is no Gauss loop. The determinant of the Jacobian supplies the volume scaling, and the factor of 2.0 comes from the library’s reference mapping for the tetrahedron formulation.

`qeT4e` uses the same geometric path, but now the displacement vector `Ue` enters the computation:

```python
Ee = (B @ Ue).reshape(-1)
Se = Ee @ D
qe = (B.T @ Se.reshape(-1, 1)) * np.linalg.det(J)
```

This is the complete local recovery path. There is no hidden state between calls; every quantity is rebuilt from the inputs.

### 1.3 Batch assembly

The global routines `kT4e` and `qT4e` are the batch versions used by the assembled solver. They do three things:

1. Extract the element node indices from `T`.
2. Gather `Xe` and the per-element material rows.
3. Scatter the local element contributions back into the global arrays.

The implementation uses `np.einsum` and `np.linalg.solve` in batch mode, so the solver can process many T4 elements without a Python loop. The scatter step uses `element_dof_indices(..., one_based=False)` plus `np.add.at`, which keeps the assembly logic explicit and avoids manual indexing bugs.

## 2. The 8-Node Hexahedron Element

The H8 element is the library’s trilinear brick. It also has 3 translational degrees of freedom per node, so one element contributes 24 local degrees of freedom. Unlike T4, the derivatives are not constant, so the element must be integrated at Gauss points.

### 2.1 Trilinear interpolation

The code uses the standard tensor-product interpolation on the parent cube. The derivative helper `_hexa_dN_batch` returns the parent-space gradients for all Gauss points at once. That is the main reason the H8 path is still vectorized even though the element is not constant strain.

The geometry enters through the nodal coordinates `Xe`. The implementation never stores a per-element state object. It recomputes everything from `Xe` and `Ge` on demand.

### 2.2 Gauss integration path

`keh8e` uses the standard 2 x 2 x 2 rule. The Gauss points are defined once:

```python
gauss_points = np.array(
    [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
    ],
    dtype=float,
) / np.sqrt(3.0)
```

Then the code evaluates the Jacobian for all points in one shot:

```python
dN = _hexa_dN_batch(gauss_points)
Jt = np.einsum("gik,kj->gij", dN, Xe)
dN_global = np.linalg.solve(Jt, dN)
B = _solid_B_batch(dN_global)
```

That sequence is the core of the element. `Jt` carries the geometry, `dN_global` carries the mapped gradients, and `B` carries the strain-displacement relation used by the constitutive law.

The stiffness matrix is then accumulated with a single Einstein summation:

```python
return np.einsum(
    "g,gik,kl,glj->ij",
    np.linalg.det(Jt),
    B.transpose(0, 2, 1),
    D,
    B,
)
```

This is the same bilinear form as the textbook integral, but expressed in the batch layout that the code actually uses.

### 2.3 Internal forces and stress recovery

`qeh8e` follows the same Gauss-point path, but also multiplies by the element displacement vector `Ue`:

```python
Ee = np.einsum("gij,jk->gi", B, Ue).reshape(8, 6)
Se = np.einsum("gi,ij->gj", Ee, D)
qe = np.einsum(
    "g,gij,gj->i",
    np.linalg.det(Jt),
    B.transpose(0, 2, 1),
    Se,
).reshape(-1, 1)
```

So the state flow is:

`Xe -> dN -> J -> dN_global -> B -> Ee -> Se -> qe`

The batch assemblers `kh8e` and `qh8e` then scatter those element results into the global arrays. The code also keeps an explicit guard against unsupported 20-node hexahedra, which matches the current repository limitation.

## 3. Practical reading guide

If you are tracing the implementation, start with `keT4e`, `qeT4e`, `keh8e`, and `qeh8e`. Those are the pure element kernels. Then read `kT4e`, `qT4e`, `kh8e`, and `qh8e` to see how the library turns local matrices into global assembly. The important design choice is that state is not cached inside the element routines; every call recomputes the geometry from the provided coordinates and material row.

