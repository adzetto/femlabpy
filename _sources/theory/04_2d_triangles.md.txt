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

# 2D Triangular Elements (CST)

This chapter covers the mathematical foundations, finite element formulation, and Python implementation of the 3-node triangular element. Often referred to as the Constant Strain Triangle (CST) in structural mechanics, it is one of the earliest and simplest finite elements developed. Despite its simplicity, it is robust and forms the foundational stepping stone for more complex isoparametric element derivations.

In this text, we will focus on both the theoretical derivation using area coordinates and the practical, highly optimized vectorized implementation available in `femlabpy`.

## Area Coordinates and Shape Functions

In a 1D line element, we use length-based natural coordinates. In a 2D triangular element, it is most natural to use **area coordinates** (also known as barycentric coordinates). 

Consider a triangle with vertices $(x_1, y_1)$, $(x_2, y_2)$, and $(x_3, y_3)$. Its total area $A$ can be computed using the determinant of the Jacobian:

$$
2A = \det \begin{bmatrix} 1 & x_1 & y_1 \\ 1 & x_2 & y_2 \\ 1 & x_3 & y_3 \end{bmatrix} = (x_2 y_3 - x_3 y_2) + (x_3 y_1 - x_1 y_3) + (x_1 y_2 - x_2 y_1)
$$

For any interior point $P(x, y)$, we can form three sub-triangles by connecting $P$ to the vertices. Let the areas of these sub-triangles be $A_1$, $A_2$, and $A_3$ (where $A_i$ is the area of the sub-triangle opposite to vertex $i$). The area coordinates $L_i$ are defined as the ratio of these sub-areas to the total area:

$$
L_1 = \frac{A_1}{A}, \quad L_2 = \frac{A_2}{A}, \quad L_3 = \frac{A_3}{A}
$$

By definition, these coordinates must satisfy the partition of unity:

$$
L_1 + L_2 + L_3 = 1
$$

In a 3-node triangular element, the shape functions $N_i(x,y)$ are exactly the area coordinates! Thus:

$$
N_1 = L_1, \quad N_2 = L_2, \quad N_3 = L_3
$$

The relationship between Cartesian coordinates $(x,y)$ and area coordinates $(L_1, L_2, L_3)$ is linear:

$$
\begin{bmatrix} 1 \\ x \\ y \end{bmatrix} = \begin{bmatrix} 1 & 1 & 1 \\ x_1 & x_2 & x_3 \\ y_1 & y_2 & y_3 \end{bmatrix} \begin{bmatrix} L_1 \\ L_2 \\ L_3 \end{bmatrix}
$$

Inverting this relationship gives the explicit polynomial form for the shape functions:

$$
N_i(x, y) = \frac{1}{2A} (a_i + b_i x + c_i y)
$$

where the geometric constants are derived from the cyclic permutations of the node indices (e.g., for $i=1$, $j=2$, $k=3$):

$$
a_i = x_j y_k - x_k y_j, \quad b_i = y_j - y_k, \quad c_i = x_k - x_j
$$

The gradients of the shape functions are therefore purely constant over the element:

$$
\frac{\partial N_i}{\partial x} = \frac{b_i}{2A}, \quad \frac{\partial N_i}{\partial y} = \frac{c_i}{2A}
$$

This fact gives the Constant Strain Triangle its name: because the displacement gradients (strains) depend on the shape function derivatives, the strain field is uniform across the entire triangle.

## Solid Mechanics: The CST Element (`ket3e` and `qet3e`)

In 2D plane elasticity (either plane stress or plane strain), each node has two degrees of freedom (DOFs): translation in $x$ and $y$. For a 3-node element, the element displacement vector is:

$$
\mathbf{u}^e = \begin{bmatrix} u_1 & v_1 & u_2 & v_2 & u_3 & v_3 \end{bmatrix}^T
$$

### The Strain-Displacement Matrix $\mathbf{B}$

The strains $\boldsymbol{\varepsilon} = [\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]^T$ are obtained by differentiating the assumed displacement field:

$$
\boldsymbol{\varepsilon} = \begin{bmatrix} \frac{\partial u}{\partial x} \\ \frac{\partial v}{\partial y} \\ \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \end{bmatrix} = \mathbf{B} \mathbf{u}^e
$$

Substituting the shape function derivatives, the exact strain-displacement matrix $\mathbf{B}$ for the T3 element becomes a constant $3 \times 6$ matrix:

$$
\mathbf{B} = \frac{1}{2A} \begin{bmatrix} 
b_1 & 0 & b_2 & 0 & b_3 & 0 \\
0 & c_1 & 0 & c_2 & 0 & c_3 \\
c_1 & b_1 & c_2 & b_2 & c_3 & b_3
\end{bmatrix}
$$

### Element Stiffness Matrix

The element stiffness matrix is derived from the principle of virtual work:

$$
\mathbf{k}^e = \int_{\Omega_e} \mathbf{B}^T \mathbf{D} \mathbf{B} t \, dA
$$

where $\mathbf{D}$ is the constitutive (elasticity) matrix for plane stress or plane strain, and $t$ is the element thickness. Because $\mathbf{B}$ and $\mathbf{D}$ are constant, the integration over the triangle area simply multiplies the integrand by the area $A$. Assuming unit thickness $t=1$:

$$
\mathbf{k}^e = \mathbf{B}^T \mathbf{D} \mathbf{B} A
$$

In `femlabpy`, this matrix is computed by the `ket3e` function. The code extracts the geometric differences efficiently using a helper function `_triangle_geometry`, constructs the $\mathbf{B}$ matrix, calls `_elastic_matrix` for the $\mathbf{D}$ tensor, and finally performs the matrix multiplication.

```python
import numpy as np
from femlabpy.elements.triangles import ket3e

# Nodal coordinates for a right triangle
Xe = np.array([[0.0, 0.0],
               [1.0, 0.0],
               [0.0, 1.0]])

# Material props: [E, nu, type] (type 1 = plane stress)
Ge = np.array([200e9, 0.3, 1])

# Generate the 6x6 element stiffness matrix
Ke = ket3e(Xe, Ge)
print(Ke.shape)
```

### Stress Recovery

Because the strain is constant, the stresses are also constant within each CST element. We evaluate them at the centroid of the element (though realistically, it represents the entire domain $\Omega_e$). The element internal forces $\mathbf{q}^e$ can be recovered by:

$$
\boldsymbol{\sigma} = \mathbf{D} \mathbf{B} \mathbf{u}^e, \quad \mathbf{q}^e = \mathbf{B}^T \boldsymbol{\sigma} A
$$

The `femlabpy` function `qet3e` performs exactly this recovery. It accepts the element displacement vector `Ue` and returns the internal forces, the local stress components, and the strain tensor.

## Potential Flow and Heat Transfer (`ket3p`)

Beyond solid mechanics, the T3 element is widely used for scalar field problems governed by the Poisson or Laplace equation, such as steady-state heat conduction, groundwater flow, or ideal fluid potential flow.

In these problems, each node has a single scalar degree of freedom (e.g., temperature $T$ or hydraulic head $\phi$). The primary unknown is a scalar field, and the flux $\mathbf{q}$ is proportional to the gradient of this field:

$$
\mathbf{q} = -k \nabla T
$$

where $k$ is the conductivity (or permeability).

The element conductivity matrix $\mathbf{k}_p^e$ is derived from the weak form:

$$
\mathbf{k}_p^e = \int_{\Omega_e} (\nabla \mathbf{N})^T k (\nabla \mathbf{N}) \, dA
$$

For the 3-node triangle, the gradient operator matrix $\mathbf{B}_p$ is a $2 \times 3$ matrix:

$$
\mathbf{B}_p = \nabla \mathbf{N} = \frac{1}{2A} \begin{bmatrix} b_1 & b_2 & b_3 \\ c_1 & c_2 & c_3 \end{bmatrix}
$$

Evaluating the integral (assuming constant $k$ and isotropic material $\mathbf{D} = k\mathbf{I}$), we get:

$$
\mathbf{k}_p^e = \mathbf{B}_p^T (k \mathbf{I}) \mathbf{B}_p A
$$

If a volumetric source term $b$ (such as internal heat generation) is present, its contribution to the element load vector involves the integral of the shape functions. For a Galerkin formulation, this generates a term added to the matrix if structured as a reaction/source boundary value, evaluated exactly over the triangle area using area coordinates:

$$
\int_{\Omega_e} N_i N_j \, dA = \frac{A}{12} (1 + \delta_{ij})
$$

In `femlabpy`, `ket3p` implements this formulation:

```python
from femlabpy.elements.triangles import ket3p

# Material properties: [conductivity, source_term]
Ge_thermal = np.array([50.0, 0.0])

# Generate the 3x3 conductivity matrix
Kp = ket3p(Xe, Ge_thermal)
print(Kp.shape)
```

Likewise, `qet3p` handles the recovery of the gradients and fluxes for these scalar elements.

## Vectorized Assembly (`kt3e` and `qt3p`)

In realistic engineering simulations, models contain tens of thousands to millions of triangular elements. Looping over every element in Python is slow due to interpreter overhead. `femlabpy` overcomes this bottleneck by vectorizing the finite element assembly.

Functions like `kt3e` (solid stiffness) and `kt3p` (potential conductivity) evaluate **all** elements simultaneously. 

### Geometric Vectorization
Instead of computing the geometry for one triangle at a time, `_triangle_batch_geometry` slices a 3D tensor of coordinates `Xe` of shape `(N, 3, 2)`. It computes the edge differences and areas for all `N` elements efficiently using `numpy` array operations:

```python
edges = np.stack([Xe[:, 2] - Xe[:, 1], Xe[:, 0] - Xe[:, 2], Xe[:, 1] - Xe[:, 0]], axis=1)
area = 0.5 * np.abs(
    (Xe[:, 1, 0] - Xe[:, 0, 0]) * (Xe[:, 2, 1] - Xe[:, 0, 1])
    - (Xe[:, 2, 0] - Xe[:, 0, 0]) * (Xe[:, 1, 1] - Xe[:, 0, 1])
)
```

### Tensor Contraction (Einsum)
Once the batched $\mathbf{B}$ tensor (shape `(N, 3, 6)`) and constitutive tensor $\mathbf{D}$ (shape `(N, 3, 3)`) are assembled, the element stiffness matrices for the entire domain are calculated via a single `np.einsum` tensor contraction:

$$
\mathbf{K}_{eij} = A_e \sum_{k,l} B_{eki}^T D_{ekl} B_{elj}
$$

In code:
```python
element_matrices = area[:, None, None] * np.einsum(
    "eik,ekl,elj->eij", B.transpose(0, 2, 1), D, B
)
```

This single command utilizes highly optimized C/Fortran BLAS backends via NumPy, bypassing the Python interpreter entirely. 

### Scatter to Global Matrix
After generating the element matrices, `femlabpy` maps the local node contributions to the global sparse stiffness matrix. The mapping uses `element_dof_indices` to generate the correct row and column arrays. If `scipy.sparse` is available, it constructs a Coordinate format (`COO`) matrix and adds it to the global matrix:

```python
delta = sp.coo_matrix(
    (element_matrices.reshape(-1), (scatter_rows, scatter_cols)),
    shape=K.shape,
    dtype=float,
)
K = (K.tocsr() + delta.tocsr()).tolil()
```

This vectorized approach results in orders-of-magnitude faster assembly compared to pure-Python `for` loops, rendering `femlabpy` exceptionally capable of executing large-scale structural and thermal analyses natively in Python.

## Summary

The Constant Strain Triangle is structurally elementary but provides the foundation for 2D finite element analysis. By exploiting exact integration with area coordinates, constant field gradients, and extreme vectorization of tensor operations via `einsum`, `femlabpy` achieves high-performance simulation standards suitable for both academic study and rigorous computational modeling.
