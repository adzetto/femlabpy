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

This chapter derives the mathematics behind the 3D solid elements in `femlabpy.elements.solids`. We begin with the fundamental strain-displacement relations in three dimensions, then specialize to the 4-node tetrahedron (T4) and 8-node hexahedron (H8).

## Preliminaries: 3D Elasticity

### Strain tensor

For a displacement field $\mathbf{u} = [u, v, w]^T$ in three dimensions, the infinitesimal strain tensor under small deformation assumptions is:

$$
\boldsymbol{\varepsilon} = \frac{1}{2}\left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \right)
$$

Expanding this in Cartesian coordinates yields six independent strain components. Using Voigt notation, we arrange them as:

$$
\boldsymbol{\varepsilon} = \begin{bmatrix} \varepsilon_{xx} \\ \varepsilon_{yy} \\ \varepsilon_{zz} \\ \gamma_{xy} \\ \gamma_{yz} \\ \gamma_{zx} \end{bmatrix} = \begin{bmatrix} \frac{\partial u}{\partial x} \\ \frac{\partial v}{\partial y} \\ \frac{\partial w}{\partial z} \\ \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \\ \frac{\partial v}{\partial z} + \frac{\partial w}{\partial y} \\ \frac{\partial w}{\partial x} + \frac{\partial u}{\partial z} \end{bmatrix}
$$

### Constitutive relation

For isotropic linear elasticity, the stress-strain relation is:

$$
\boldsymbol{\sigma} = \mathbf{D} \boldsymbol{\varepsilon}
$$

where the $6 \times 6$ elasticity matrix $\mathbf{D}$ is:

$$
\mathbf{D} = \frac{E}{(1+\nu)(1-2\nu)} \begin{bmatrix} 1-\nu & \nu & \nu & 0 & 0 & 0 \\ \nu & 1-\nu & \nu & 0 & 0 & 0 \\ \nu & \nu & 1-\nu & 0 & 0 & 0 \\ 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 & 0 \\ 0 & 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 \\ 0 & 0 & 0 & 0 & 0 & \frac{1-2\nu}{2} \end{bmatrix}
$$

This matrix is implemented in `_elastic3d_matrix(Ge)` where `Ge = [E, nu]`.

---

## 1. The 4-Node Tetrahedron Element (T4)

```{figure} ../_static/figures/fig_t4_element.png
:align: center
:width: 480px

**Figure 6.1:** 4-node tetrahedron (T4) element with 3 translational DOFs per node. Dashed edge indicates hidden line.
```

The T4 element has 4 nodes with 3 DOFs each, giving 12 element DOFs total. The topology row is `[n1, n2, n3, n4, mat_id]`.

### 1.1 Shape functions

The T4 uses volume (barycentric) coordinates $L_1, L_2, L_3, L_4$ as shape functions. For a point $\mathbf{x}$ inside the tetrahedron with vertices $\mathbf{x}_1, \mathbf{x}_2, \mathbf{x}_3, \mathbf{x}_4$:

$$
N_i(\mathbf{x}) = L_i = \frac{V_i}{V}
$$

where $V$ is the total volume and $V_i$ is the volume of the sub-tetrahedron formed by replacing vertex $i$ with point $\mathbf{x}$.

The volume of a tetrahedron with vertices at $\mathbf{x}_1, \mathbf{x}_2, \mathbf{x}_3, \mathbf{x}_4$ is:

$$
V = \frac{1}{6} \det \begin{bmatrix} x_2 - x_1 & x_3 - x_1 & x_4 - x_1 \\ y_2 - y_1 & y_3 - y_1 & y_4 - y_1 \\ z_2 - z_1 & z_3 - z_1 & z_4 - z_1 \end{bmatrix} = \frac{1}{6} \det(\mathbf{J})
$$

Since $L_1 + L_2 + L_3 + L_4 = 1$, we can parameterize with three independent coordinates. The parent-space derivatives are constant:

$$
\frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}} = \begin{bmatrix} \frac{\partial N_1}{\partial \xi} & \frac{\partial N_2}{\partial \xi} & \frac{\partial N_3}{\partial \xi} & \frac{\partial N_4}{\partial \xi} \\ \frac{\partial N_1}{\partial \eta} & \frac{\partial N_2}{\partial \eta} & \frac{\partial N_3}{\partial \eta} & \frac{\partial N_4}{\partial \eta} \\ \frac{\partial N_1}{\partial \zeta} & \frac{\partial N_2}{\partial \zeta} & \frac{\partial N_3}{\partial \zeta} & \frac{\partial N_4}{\partial \zeta} \end{bmatrix} = \begin{bmatrix} 1 & 0 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & -1 \end{bmatrix}
$$

### 1.2 Jacobian and physical derivatives

The Jacobian matrix maps parent coordinates to physical coordinates:

$$
\mathbf{J} = \frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}} \cdot \mathbf{X}_e
$$

where $\mathbf{X}_e$ is the $4 \times 3$ nodal coordinate matrix. The physical-space shape function derivatives are:

$$
\frac{\partial \mathbf{N}}{\partial \mathbf{x}} = \mathbf{J}^{-1} \frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}}
$$

In code, this is computed via `np.linalg.solve(J, dN)` to avoid explicit inversion.

### 1.3 Strain-displacement matrix

For 3D solids, the element displacement vector is:

$$
\mathbf{u}^e = [u_1, v_1, w_1, u_2, v_2, w_2, u_3, v_3, w_3, u_4, v_4, w_4]^T
$$

The $6 \times 12$ strain-displacement matrix $\mathbf{B}$ relates strain to nodal displacements:

$$
\boldsymbol{\varepsilon} = \mathbf{B} \mathbf{u}^e
$$

For node $i$ with shape function derivatives $\frac{\partial N_i}{\partial x}, \frac{\partial N_i}{\partial y}, \frac{\partial N_i}{\partial z}$, the contribution to $\mathbf{B}$ is:

$$
\mathbf{B}_i = \begin{bmatrix} \frac{\partial N_i}{\partial x} & 0 & 0 \\ 0 & \frac{\partial N_i}{\partial y} & 0 \\ 0 & 0 & \frac{\partial N_i}{\partial z} \\ \frac{\partial N_i}{\partial y} & \frac{\partial N_i}{\partial x} & 0 \\ 0 & \frac{\partial N_i}{\partial z} & \frac{\partial N_i}{\partial y} \\ \frac{\partial N_i}{\partial z} & 0 & \frac{\partial N_i}{\partial x} \end{bmatrix}
$$

The full matrix is $\mathbf{B} = [\mathbf{B}_1 | \mathbf{B}_2 | \mathbf{B}_3 | \mathbf{B}_4]$.

### 1.4 Element stiffness matrix

The element stiffness matrix is derived from the principle of virtual work. The strain energy is:

$$
U = \frac{1}{2} \int_V \boldsymbol{\varepsilon}^T \boldsymbol{\sigma} \, dV = \frac{1}{2} \int_V \boldsymbol{\varepsilon}^T \mathbf{D} \boldsymbol{\varepsilon} \, dV
$$

Substituting $\boldsymbol{\varepsilon} = \mathbf{B} \mathbf{u}^e$:

$$
U = \frac{1}{2} (\mathbf{u}^e)^T \underbrace{\int_V \mathbf{B}^T \mathbf{D} \mathbf{B} \, dV}_{\mathbf{K}^e} \mathbf{u}^e
$$

Since $\mathbf{B}$ and $\mathbf{D}$ are constant over the T4 element, the integral simplifies:

$$
\mathbf{K}^e = \mathbf{B}^T \mathbf{D} \mathbf{B} \cdot V = \mathbf{B}^T \mathbf{D} \mathbf{B} \cdot \frac{1}{6} |\det(\mathbf{J})|
$$

In code, the factor appears as `2.0 * det(J)` due to the reference tetrahedron mapping convention.

### 1.5 Internal force vector

For a given displacement $\mathbf{u}^e$, the internal force vector is:

$$
\mathbf{q}^e = \int_V \mathbf{B}^T \boldsymbol{\sigma} \, dV = \int_V \mathbf{B}^T \mathbf{D} \boldsymbol{\varepsilon} \, dV
$$

With constant strain:

$$
\mathbf{q}^e = \mathbf{B}^T \mathbf{D} (\mathbf{B} \mathbf{u}^e) \cdot V
$$

---

## 2. The 8-Node Hexahedron Element (H8)

The H8 element is the trilinear brick with 8 nodes and 24 DOFs.Unlike T4, the strain field varies within the element, requiring numerical integration.

### 2.1 Trilinear shape functions

The parent domain is the cube $[-1,1]^3$ with natural coordinates $(\xi, \eta, \zeta)$. The eight shape functions are tensor products of 1D linear functions:

$$
N_i(\xi, \eta, \zeta) = \frac{1}{8}(1 + \xi_i \xi)(1 + \eta_i \eta)(1 + \zeta_i \zeta)
$$

where $(\xi_i, \eta_i, \zeta_i)$ are the corner coordinates of node $i$:

| Node | $\xi_i$ | $\eta_i$ | $\zeta_i$ |
|------|---------|----------|-----------|
| 1    | -1      | -1       | -1        |
| 2    | +1      | -1       | -1        |
| 3    | +1      | +1       | -1        |
| 4    | -1      | +1       | -1        |
| 5    | -1      | -1       | +1        |
| 6    | +1      | -1       | +1        |
| 7    | +1      | +1       | +1        |
| 8    | -1      | +1       | +1        |

### 2.2 Shape function derivatives

The parent-space derivatives are:

$$
\frac{\partial N_i}{\partial \xi} = \frac{\xi_i}{8}(1 + \eta_i \eta)(1 + \zeta_i \zeta)
$$

$$
\frac{\partial N_i}{\partial \eta} = \frac{\eta_i}{8}(1 + \xi_i \xi)(1 + \zeta_i \zeta)
$$

$$
\frac{\partial N_i}{\partial \zeta} = \frac{\zeta_i}{8}(1 + \xi_i \xi)(1 + \eta_i \eta)
$$

These vary with position, unlike T4.

### 2.3 Isoparametric mapping

Physical coordinates are interpolated from nodal coordinates:

$$
\mathbf{x}(\xi, \eta, \zeta) = \sum_{i=1}^{8} N_i(\xi, \eta, \zeta) \mathbf{x}_i
$$

The Jacobian matrix at each point is:

$$
\mathbf{J} = \begin{bmatrix} \frac{\partial x}{\partial \xi} & \frac{\partial y}{\partial \xi} & \frac{\partial z}{\partial \xi} \\ \frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta} & \frac{\partial z}{\partial \eta} \\ \frac{\partial x}{\partial \zeta} & \frac{\partial y}{\partial \zeta} & \frac{\partial z}{\partial \zeta} \end{bmatrix} = \frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}} \cdot \mathbf{X}_e
$$

Physical derivatives are obtained by solving:

$$
\frac{\partial \mathbf{N}}{\partial \mathbf{x}} = \mathbf{J}^{-1} \frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}}
$$

### 2.4 Gauss-Legendre integration

Since the integrand varies within the element, we use numerical quadrature:

$$
\mathbf{K}^e = \int_{-1}^{1} \int_{-1}^{1} \int_{-1}^{1} \mathbf{B}^T \mathbf{D} \mathbf{B} \, |\det(\mathbf{J})| \, d\xi \, d\eta \, d\zeta
$$

The $2 \times 2 \times 2$ Gauss rule uses 8 points at $\xi, \eta, \zeta = \pm 1/\sqrt{3}$ with unit weights:

$$
\mathbf{K}^e = \sum_{g=1}^{8} w_g \, \mathbf{B}_g^T \mathbf{D} \mathbf{B}_g \, |\det(\mathbf{J}_g)|
$$

For the standard 2-point rule, $w_g = 1$ for all points.

### 2.5 Stiffness matrix computation

At each Gauss point $g$:

1. Evaluate parent derivatives $\frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}}\big|_g$
2. Compute Jacobian $\mathbf{J}_g = \frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}}\big|_g \cdot \mathbf{X}_e$
3. Solve for physical derivatives $\frac{\partial \mathbf{N}}{\partial \mathbf{x}}\big|_g = \mathbf{J}_g^{-1} \frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}}\big|_g$
4. Build $\mathbf{B}_g$ from the physical derivatives
5. Accumulate $\mathbf{K}^e \mathrel{+}= \mathbf{B}_g^T \mathbf{D} \mathbf{B}_g \cdot |\det(\mathbf{J}_g)|$

In code, this is vectorized using `np.einsum` over all Gauss points simultaneously.

### 2.6 Stress recovery

At each Gauss point, strain and stress are:

$$
\boldsymbol{\varepsilon}_g = \mathbf{B}_g \mathbf{u}^e, \quad \boldsymbol{\sigma}_g = \mathbf{D} \boldsymbol{\varepsilon}_g
$$

The internal force vector is assembled as:

$$
\mathbf{q}^e = \sum_{g=1}^{8} \mathbf{B}_g^T \boldsymbol{\sigma}_g \, |\det(\mathbf{J}_g)|
$$

---

## 3. Implementation Notes

### Code structure

| Function | Description |
|----------|-------------|
| `keT4e` | T4 element stiffness (single element) |
| `qeT4e` | T4 internal force and stress recovery |
| `kT4e` | T4 batch assembly to global stiffness |
| `qT4e` | T4 batch internal force assembly |
| `keh8e` | H8 element stiffness (single element) |
| `qeh8e` | H8 internal force and stress recovery |
| `kh8e` | H8 batch assembly to global stiffness |
| `qh8e` | H8 batch internal force assembly |

### Design principles

1. **No cached state**: Every call recomputes geometry from `Xe` and material from `Ge`
2. **Vectorized loops**: Batch routines use `np.einsum` instead of Python loops
3. **Explicit indexing**: DOF scatter uses `element_dof_indices` plus `np.add.at`
