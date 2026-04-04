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

In computational mechanics and finite element analysis (FEA), three-dimensional solid elements are essential for modeling general volumetric bodies where simplifications like plane stress, plane strain, or axisymmetric assumptions are invalid. 

This chapter provides a comprehensive overview of the theoretical foundations and implementation details for two workhorse 3D elements: the 4-node linear tetrahedron (T4) and the 8-node trilinear hexahedron (H8). We will cover their shape functions, formulation of the Jacobian, construction of the strain-displacement matrix ($\mathbf{B}$), stiffness matrix ($\mathbf{K}^e$), and internal force vectors.

## 1. The 4-Node Tetrahedron Element (T4)

The 4-node tetrahedron (often called the constant strain tetrahedron, or CST in 3D) is the simplest 3D solid element. It has four nodes, each with three translational degrees of freedom ($u, v, w$), yielding a total of 12 degrees of freedom per element.

## 1.1 Volume Coordinates and Shape Functions

For a tetrahedron, it is highly convenient to define shape functions in terms of volume coordinates, $L_1, L_2, L_3, L_4$. A point $P$ inside the tetrahedron partitions the total volume $V$ into four sub-volumes $V_1, V_2, V_3, V_4$ subtended by the faces opposite to nodes 1, 2, 3, and 4. The volume coordinates are defined as:

$$
L_i = \frac{V_i}{V}, \quad i = 1, 2, 3, 4
$$

By definition, these coordinates sum to unity: $L_1 + L_2 + L_3 + L_4 = 1$.

In an isoparametric mapping, we can choose three independent natural coordinates $\xi, \eta, \zeta$ such that:
$$
\xi = L_1, \quad \eta = L_2, \quad \zeta = L_3 \implies L_4 = 1 - \xi - \eta - \zeta
$$

The shape functions $N_i$ for the T4 element are simply the volume coordinates:
$$
\begin{align*}
N_1(\xi, \eta, \zeta) &= \xi \\
N_2(\xi, \eta, \zeta) &= \eta \\
N_3(\xi, \eta, \zeta) &= \zeta \\
N_4(\xi, \eta, \zeta) &= 1 - \xi - \eta - \zeta
\end{align*}
$$

The local derivatives of the shape functions with respect to the natural coordinates are easily evaluated:
$$
\frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}} =
\begin{bmatrix}
\frac{\partial N_1}{\partial \xi} & \frac{\partial N_2}{\partial \xi} & \frac{\partial N_3}{\partial \xi} & \frac{\partial N_4}{\partial \xi} \\
\frac{\partial N_1}{\partial \eta} & \frac{\partial N_2}{\partial \eta} & \frac{\partial N_3}{\partial \eta} & \frac{\partial N_4}{\partial \eta} \\
\frac{\partial N_1}{\partial \zeta} & \frac{\partial N_2}{\partial \zeta} & \frac{\partial N_3}{\partial \zeta} & \frac{\partial N_4}{\partial \zeta}
\end{bmatrix}
=
\begin{bmatrix}
1 & 0 & 0 & -1 \\
0 & 1 & 0 & -1 \\
0 & 0 & 1 & -1
\end{bmatrix}
$$

## 1.2 The Jacobian and the B-Matrix

The Jacobian matrix $\mathbf{J}$ maps the natural coordinates to the physical Cartesian coordinates $(x, y, z)$. It is computed as:
$$
\mathbf{J} = \frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}} \mathbf{X}_e
$$
where $\mathbf{X}_e$ is the $4 \times 3$ matrix of nodal coordinates. Because the derivatives of the shape functions are constant, $\mathbf{J}$ is constant over the element. 

Using the inverse of the Jacobian, we map the shape function derivatives from natural to global coordinates:
$$
\begin{bmatrix} \frac{\partial N_i}{\partial x} \\ \frac{\partial N_i}{\partial y} \\ \frac{\partial N_i}{\partial z} \end{bmatrix} = \mathbf{J}^{-1} \begin{bmatrix} \frac{\partial N_i}{\partial \xi} \\ \frac{\partial N_i}{\partial \eta} \\ \frac{\partial N_i}{\partial \zeta} \end{bmatrix}
$$

The strain-displacement matrix $\mathbf{B}$ relates the nodal displacements $\mathbf{u}_e$ to the 6 strain components $\boldsymbol{\epsilon} = [\epsilon_{xx}, \epsilon_{yy}, \epsilon_{zz}, \gamma_{xy}, \gamma_{yz}, \gamma_{zx}]^T$. The $\mathbf{B}$ matrix is of size $6 \times 12$ and can be written as:
$$
\mathbf{B} = \begin{bmatrix} \mathbf{B}_1 & \mathbf{B}_2 & \mathbf{B}_3 & \mathbf{B}_4 \end{bmatrix}
$$
where for each node $i$:
$$
\mathbf{B}_i = \begin{bmatrix}
N_{i,x} & 0 & 0 \\
0 & N_{i,y} & 0 \\
0 & 0 & N_{i,z} \\
N_{i,y} & N_{i,x} & 0 \\
0 & N_{i,z} & N_{i,y} \\
N_{i,z} & 0 & N_{i,x}
\end{bmatrix}
$$

## 1.3 Stiffness Matrix Formulation

For linear elasticity, the element stiffness matrix $\mathbf{K}^e$ is given by:
$$
\mathbf{K}^e = \int_{V} \mathbf{B}^T \mathbf{D} \mathbf{B} \, dV
$$
where $\mathbf{D}$ is the $6 \times 6$ constitutive elasticity matrix for an isotropic material.

Since $\mathbf{B}$ and $\mathbf{D}$ are constant for the T4 element, the integral evaluates simply to:
$$
\mathbf{K}^e = \mathbf{B}^T \mathbf{D} \mathbf{B} \int_V dV = \mathbf{B}^T \mathbf{D} \mathbf{B} \, V_e
$$
The volume of the tetrahedron $V_e$ is related to the determinant of the Jacobian matrix by $V_e = \frac{1}{6} \det(\mathbf{J})$.

## 2. The 8-Node Hexahedron Element (H8)

The 8-node hexahedron is a widespread 3D solid element that utilizes a trilinear shape function interpolation. With 3 degrees of freedom per node, the element has a total of 24 degrees of freedom.

## 2.1 Trilinear Shape Functions

The element is formulated in a reference parent coordinate system defined by $\xi, \eta, \zeta \in [-1, 1]$. The eight nodes are located at the corners of this reference cube: $(\pm 1, \pm 1, \pm 1)$.

The shape functions $N_i$ for node $i$ (with parent coordinates $\xi_i, \eta_i, \zeta_i$) are given by:
$$
N_i(\xi, \eta, \zeta) = \frac{1}{8} (1 + \xi_i \xi)(1 + \eta_i \eta)(1 + \zeta_i \zeta) \quad \text{for } i = 1, \dots, 8
$$

These trilinear shape functions are $C^0$ continuous across element boundaries, ensuring displacement compatibility.

## 2.2 Jacobian and Kinematics

Unlike the T4 element, the derivatives of the H8 shape functions are not constant; they depend linearly on the local coordinates. For example, the derivative of $N_1$ (located at $(-1, -1, -1)$) with respect to $\xi$ is:
$$
\frac{\partial N_1}{\partial \xi} = -\frac{1}{8} (1 - \eta)(1 - \zeta)
$$

The $3 \times 3$ Jacobian matrix $\mathbf{J}$ is a function of the local coordinates $\xi, \eta, \zeta$ and must be evaluated at specific points during numerical integration:
$$
\mathbf{J}(\xi, \eta, \zeta) = \sum_{i=1}^8 \begin{bmatrix} \frac{\partial N_i}{\partial \xi} \\ \frac{\partial N_i}{\partial \eta} \\ \frac{\partial N_i}{\partial \zeta} \end{bmatrix} \begin{bmatrix} x_i & y_i & z_i \end{bmatrix}
$$

The $6 \times 24$ strain-displacement matrix $\mathbf{B}$ is similarly evaluated pointwise by mapping the parent derivatives using the inverse of $\mathbf{J}(\xi, \eta, \zeta)$.

## 2.3 Gauss Quadrature Integration

The element stiffness matrix is:
$$
\mathbf{K}^e = \int_{-1}^1 \int_{-1}^1 \int_{-1}^1 \mathbf{B}^T \mathbf{D} \mathbf{B} \det(\mathbf{J}) \, d\xi \, d\eta \, d\zeta
$$
Since $\mathbf{B}$ and $\mathbf{J}$ are functions of $\xi, \eta, \zeta$, we use Gauss-Legendre quadrature to numerically evaluate the integral. The standard rule for an H8 element is a $2 \times 2 \times 2$ integration scheme, which exactly integrates polynomials up to degree 3. 

The eight Gauss points are situated at $\pm 1/\sqrt{3}$ along each local axis:
$$
r_k \in \left\{-\frac{1}{\sqrt{3}}, \frac{1}{\sqrt{3}}\right\}
$$
Each integration point has a weight $W = 1 \times 1 \times 1 = 1$. The numerical stiffness matrix is:
$$
\mathbf{K}^e = \sum_{g=1}^8 \mathbf{B}_g^T \mathbf{D} \mathbf{B}_g \det(\mathbf{J}_g)
$$

## 3. Implementation in Python

In this section, we examine the practical implementation of the theoretical concepts above using Python. These functions demonstrate computation of the stiffness matrices and recovery of internal forces.

## 3.1 `keT4e` and `qeT4e`

The `keT4e` function generates the $12 \times 12$ local element stiffness matrix for the T4 element. Notice the fixed `dN` array which matches our theoretical constant gradient matrix.

```python
import numpy as np

def keT4e(Xe, Ge):
    # Xe: Nodal coordinates (4, 3)
    # Ge: Material parameters [E, nu]
    dN = np.array(
        [[1.0, 0.0, 0.0, -1.0], 
         [0.0, 1.0, 0.0, -1.0], 
         [0.0, 0.0, 1.0, -1.0]], dtype=float
    )
    J = dN @ Xe
    dN_global = np.linalg.solve(J, dN)
    
    B = _solid_B(dN_global)
    D = _elastic3d_matrix(Ge)
    
    # K = B^T * D * B * (Volume factor via det J)
    return 2.0 * (B.T @ D @ B) * np.linalg.det(J)
```

The `qeT4e` function computes the element's internal force vector $\mathbf{q}_e = \mathbf{K}^e \mathbf{u}_e$, as well as stresses ($\boldsymbol{\sigma}$) and strains ($\boldsymbol{\epsilon}$).

```python
def qeT4e(Xe, Ge, Ue):
    dN = np.array(
        [[1.0, 0.0, 0.0, -1.0], [0.0, 1.0, 0.0, -1.0], [0.0, 0.0, 1.0, -1.0]],
        dtype=float,
    )
    J = dN @ Xe
    dN_global = np.linalg.solve(J, dN)
    
    B = _solid_B(dN_global)
    D = _elastic3d_matrix(Ge)
    
    Ee = (B @ Ue).reshape(-1)
    Se = Ee @ D
    qe = (B.T @ Se.reshape(-1, 1)) * np.linalg.det(J)
    
    return qe, Se, Ee
```

## 3.2 `keh8e` and `qeh8e`

The `keh8e` implementation performs numerical integration over the 8 Gauss points to assemble the $24 \times 24$ stiffness matrix. Vectorization (`_hexa_dN_batch`) is often used to rapidly evaluate quantities over all Gauss points simultaneously.

```python
def keh8e(Xe, Ge):
    D = _elastic3d_matrix(Ge)
    
    # 2x2x2 Gauss Quadrature Points
    gauss_points = np.array([
        [-1., -1., -1.], [1., -1., -1.], [1., 1., -1.], [-1., 1., -1.],
        [-1., -1., 1.],  [1., -1., 1.],  [1., 1., 1.],  [-1., 1., 1.]
    ], dtype=float) / np.sqrt(3.0)
    
    # Evaluate derivatives of N at all 8 Gauss points
    dN = _hexa_dN_batch(gauss_points)
    
    # Compute Jacobians and inverse transformations
    Jt = np.einsum("gik,kj->gij", dN, Xe)
    dN_global = np.linalg.solve(Jt, dN)
    
    B = _solid_B_batch(dN_global)
    
    # Integral over the element domain using Einsum for contraction
    return np.einsum(
        "g,gik,kl,glj->ij",
        np.linalg.det(Jt),
        B.transpose(0, 2, 1),
        D,
        B,
    )
```

The recovery function `qeh8e` applies identical integration structures to calculate stresses and strains at each of the eight internal Gauss points:

```python
def qeh8e(Xe, Ge, Ue):
    # Setup D, Gauss points, and B matrix same as in keh8e ...
    
    # Element Strains at Gauss Points
    Ee = np.einsum("gij,jk->gi", B, Ue).reshape(8, 6)
    
    # Element Stresses at Gauss Points
    Se = np.einsum("gi,ij->gj", Ee, D)
    
    # Internal Force Vector
    qe = np.einsum(
        "g,gij,gj->i",
        np.linalg.det(Jt),
        B.transpose(0, 2, 1),
        Se,
    ).reshape(-1, 1)
    
    return qe, Se, Ee
```

## Summary

In 3D solid mechanics, tetrahedral elements (T4) offer unmatched meshing flexibility, particularly for highly complex CAD geometries. However, due to their constant strain formulation, they suffer from artificial stiffening (shear locking) and demand a remarkably fine mesh for accurate bending representations.

Hexahedral elements (H8), conversely, exhibit superior accuracy with fewer degrees of freedom due to their trilinear interpolation. The primary challenge with H8 elements rests on the difficulty of automatically meshing arbitrary 3D volumes into purely hexahedral grids. Modern solvers leverage both element types, adapting integration techniques and topologies to meet the stringent demands of computational precision.
