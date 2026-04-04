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

# Chapter 5: 2D Isoparametric Quadrilateral Elements

The finite element method (FEM) heavily relies on the concept of isoparametric elements. These elements use the same shape functions to interpolate both the geometric coordinates and the field variables (such as displacements or temperatures) within the element. This powerful concept allows elements to take on distorted, non-rectangular shapes in the global Cartesian coordinate system while being mapped from a simple, standard square shape in a local, "parent" coordinate system.

In this comprehensive chapter, we delve into the theory, mathematics, and computational implementation of the 4-node bilinear isoparametric quadrilateral element (Q4 element). We will explore the formulation of its shape functions, the Jacobian matrix used for coordinate transformation, the strain-displacement matrix, and the application of Gauss-Legendre numerical integration. Finally, we will dissect its software implementation within `femlabpy`, covering both solid mechanics (elasticity) and scalar potential flow (e.g., heat conduction, seepage).

## 1. Bilinear Shape Functions

The Q4 element has four nodes, typically numbered counter-clockwise starting from the bottom-left corner. To facilitate integration and generic element formulation, we define a "parent" element in a local coordinate system with axes denoted by $\xi$ and $\eta$. The parent element is a square with its domain defined by $\xi \in [-1, 1]$ and $\eta \in [-1, 1]$. 

The nodes of the parent element are located at:
* Node 1: $(-1, -1)$
* Node 2: $(1, -1)$
* Node 3: $(1, 1)$
* Node 4: $(-1, 1)$

For an interpolation function to be valid, it must take the value of $1$ at its corresponding node and $0$ at all other nodes. For a 4-node element, the simplest polynomials that satisfy these conditions and represent a complete bi-linear field are the bilinear shape functions. These are constructed by taking the tensor product of 1D linear Lagrange polynomials.

The general expression for the shape function $N_i$ associated with node $i$ (having coordinates $\xi_i, \eta_i$) is given by:

$$
N_i(\xi, \eta) = \frac{1}{4} (1 + \xi_i \xi) (1 + \eta_i \eta) \quad \text{for } i = 1, 2, 3, 4
$$

Expanding this for each of the four nodes, we obtain the individual shape functions:

$$
\begin{aligned}
N_1(\xi, \eta) &= \frac{1}{4} (1 - \xi) (1 - \eta) \\
N_2(\xi, \eta) &= \frac{1}{4} (1 + \xi) (1 - \eta) \\
N_3(\xi, \eta) &= \frac{1}{4} (1 + \xi) (1 + \eta) \\
N_4(\xi, \eta) &= \frac{1}{4} (1 - \xi) (1 + \eta)
\end{aligned}
$$

These shape functions possess several important properties:
1.  **Kronecker delta property:** $N_i(\xi_j, \eta_j) = \delta_{ij}$.
2.  **Partition of unity:** $\sum_{i=1}^4 N_i(\xi, \eta) = 1$ for any point $(\xi, \eta)$ within the element. This ensures that a rigid body translation is exactly represented.
3.  **Linear representation:** They can exactly represent a linear field $C_1 + C_2 x + C_3 y$.

## 2. Isoparametric Mapping and the Jacobian Matrix

In the isoparametric formulation, the global Cartesian coordinates $(x, y)$ of any point within the element are interpolated from the nodal coordinates $(x_i, y_i)$ using the same shape functions $N_i$:

$$
x(\xi, \eta) = \sum_{i=1}^4 N_i(\xi, \eta) x_i, \quad y(\xi, \eta) = \sum_{i=1}^4 N_i(\xi, \eta) y_i
$$

Similarly, the primary field variable—let's consider the displacement vector $\mathbf{u} = [u, v]^T$ for solid mechanics—is interpolated as:

$$
u(\xi, \eta) = \sum_{i=1}^4 N_i(\xi, \eta) u_i, \quad v(\xi, \eta) = \sum_{i=1}^4 N_i(\xi, \eta) v_i
$$

To compute strains or fluxes, we need the derivatives of the shape functions with respect to the global coordinates $x$ and $y$. However, the shape functions are defined in terms of the local coordinates $\xi$ and $\eta$. We must employ the chain rule of differentiation:

$$
\begin{aligned}
\frac{\partial N_i}{\partial \xi} &= \frac{\partial N_i}{\partial x} \frac{\partial x}{\partial \xi} + \frac{\partial N_i}{\partial y} \frac{\partial y}{\partial \xi} \\
\frac{\partial N_i}{\partial \eta} &= \frac{\partial N_i}{\partial x} \frac{\partial x}{\partial \eta} + \frac{\partial N_i}{\partial y} \frac{\partial y}{\partial \eta}
\end{aligned}
$$

This relationship can be expressed elegantly in matrix form using the **Jacobian matrix**, $\mathbf{J}$:

$$
\begin{Bmatrix}
\frac{\partial N_i}{\partial \xi} \\
\frac{\partial N_i}{\partial \eta}
\end{Bmatrix}
=
\begin{bmatrix}
\frac{\partial x}{\partial \xi} & \frac{\partial y}{\partial \xi} \\
\frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta}
\end{bmatrix}
\begin{Bmatrix}
\frac{\partial N_i}{\partial x} \\
\frac{\partial N_i}{\partial y}
\end{Bmatrix}
= \mathbf{J}
\begin{Bmatrix}
\frac{\partial N_i}{\partial x} \\
\frac{\partial N_i}{\partial y}
\end{Bmatrix}
$$

The components of the Jacobian matrix $\mathbf{J}$ are evaluated using the coordinate interpolation equations:

$$
\mathbf{J} = \begin{bmatrix}
\sum_{i=1}^4 \frac{\partial N_i}{\partial \xi} x_i & \sum_{i=1}^4 \frac{\partial N_i}{\partial \xi} y_i \\
\sum_{i=1}^4 \frac{\partial N_i}{\partial \eta} x_i & \sum_{i=1}^4 \frac{\partial N_i}{\partial \eta} y_i
\end{bmatrix}
$$

To find the desired Cartesian derivatives of the shape functions, we simply invert the Jacobian matrix:

$$
\begin{Bmatrix}
\frac{\partial N_i}{\partial x} \\
\frac{\partial N_i}{\partial y}
\end{Bmatrix}
= \mathbf{J}^{-1}
\begin{Bmatrix}
\frac{\partial N_i}{\partial \xi} \\
\frac{\partial N_i}{\partial \eta}
\end{Bmatrix}
$$

The determinant of the Jacobian, $|\mathbf{J}|$, provides the scale factor relating the infinitesimal area in the parent domain to the physical domain: $dx dy = |\mathbf{J}| d\xi d\eta$. For the mapping to be valid and invertible, $|\mathbf{J}|$ must be strictly positive everywhere inside the element, implying the element cannot be excessively distorted, self-intersecting, or have interior angles $\geq 180^\circ$.

## 3. The Strain-Displacement Matrix ($\mathbf{B}$)

In 2D plane elasticity (plane stress or plane strain), the strain vector is given by $\boldsymbol{\varepsilon} = [\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]^T$. Using the small-strain kinematic relations:

$$
\boldsymbol{\varepsilon} = \begin{Bmatrix}
\frac{\partial u}{\partial x} \\
\frac{\partial v}{\partial y} \\
\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}
\end{Bmatrix}
$$

Substituting the interpolated displacements into the strain definitions yields the relationship between the strains at any point and the nodal displacement vector $\mathbf{d}_e = [u_1, v_1, u_2, v_2, u_3, v_3, u_4, v_4]^T$:

$$
\boldsymbol{\varepsilon} = \mathbf{B} \mathbf{d}_e
$$

Where $\mathbf{B}$ is the strain-displacement matrix. For the Q4 element, $\mathbf{B}$ is a $3 \times 8$ matrix constructed from sub-matrices $\mathbf{B}_i$ for each node $i$:

$$
\mathbf{B} = \begin{bmatrix} \mathbf{B}_1 & \mathbf{B}_2 & \mathbf{B}_3 & \mathbf{B}_4 \end{bmatrix}
$$

Each nodal sub-matrix $\mathbf{B}_i$ is defined using the Cartesian derivatives of the shape functions:

$$
\mathbf{B}_i = \begin{bmatrix}
\frac{\partial N_i}{\partial x} & 0 \\
0 & \frac{\partial N_i}{\partial y} \\
\frac{\partial N_i}{\partial y} & \frac{\partial N_i}{\partial x}
\end{bmatrix}
$$

## 4. Gauss-Legendre Numerical Integration

The element stiffness matrix $\mathbf{K}_e$ is defined by the integral of the strain energy over the element volume $V_e$:

$$
\mathbf{K}_e = \int_{V_e} \mathbf{B}^T \mathbf{D} \mathbf{B} \, dV
$$

Assuming a uniform thickness $t$, this becomes an area integral over the physical domain $\Omega_e$:

$$
\mathbf{K}_e = t \int_{\Omega_e} \mathbf{B}^T \mathbf{D} \mathbf{B} \, dx dy
$$

Where $\mathbf{D}$ is the constitutive (material) matrix. Because the entries of $\mathbf{B}$ are rational functions (involving the inverse of the Jacobian determinant), exact analytical integration is generally impossible for an arbitrary quadrilateral. We must transform the integral to the parent domain and use numerical quadrature:

$$
\mathbf{K}_e = t \int_{-1}^{1} \int_{-1}^{1} \mathbf{B}^T(\xi, \eta) \mathbf{D} \mathbf{B}(\xi, \eta) |\mathbf{J}(\xi, \eta)| \, d\xi d\eta
$$

The standard approach is the $2 \times 2$ Gauss-Legendre quadrature. This rule uses four integration points (Gauss points) strategically located within the parent element. The integral is approximated as a weighted sum of the integrand evaluated at these points:

$$
\mathbf{K}_e \approx \sum_{i=1}^2 \sum_{j=1}^2 w_i w_j \left[ \mathbf{B}^T(\xi_i, \eta_j) \mathbf{D} \mathbf{B}(\xi_i, \eta_j) |\mathbf{J}(\xi_i, \eta_j)| t \right]
$$

For the $2 \times 2$ rule, the coordinates and weights are:
*   Points $\xi_i, \eta_j \in \left\{ -\frac{1}{\sqrt{3}}, \frac{1}{\sqrt{3}} \right\}$
*   Weights $w_i = w_j = 1.0$

This $2 \times 2$ rule exactly integrates polynomials up to degree 3, which is sufficient to accurately evaluate the stiffness matrix and prevent hourglassing (spurious zero-energy modes) that occurs if reduced $1 \times 1$ integration is used.

## 5. Software Implementation in `femlabpy`

The `femlabpy` package provides robust and efficient implementations of Q4 elements for both solid mechanics and potential flow problems. Let's examine the core functions located in `src/femlabpy/elements/quads.py`.

### 5.1 Elastic Stiffness Matrix: `keq4e`

The `keq4e` function calculates the $8 \times 8$ element stiffness matrix for a linearly elastic plane stress or plane strain problem.

```python
import numpy as np
from femlabpy.elements.quads import keq4e

# Define the nodal coordinates [x, y] for the 4 nodes
Xe = np.array([
    [0.0, 0.0],
    [2.0, 0.0],
    [2.0, 1.5],
    [0.0, 1.5]
])

# Define material properties [E, nu, type, thickness]
# type 1 = plane stress (default), type 2 = plane strain
Ge = np.array([200e9, 0.3, 1, 0.05]) # Plane stress, thickness = 0.05

# Compute the stiffness matrix
Ke = keq4e(Xe, Ge)
print("Element Stiffness Matrix Shape:", Ke.shape)
```

Internally, `keq4e` loops over the $2 \times 2$ Gauss points. At each point, it computes the local shape function derivatives (`_q4_dN`), evaluates the Jacobian (`Jt`), inverts it to get global derivatives, forms the $\mathbf{B}$ matrix (`_q4_B`), and accumulates the integrand into `Ke`.

### 5.2 Internal Forces and Stresses: `qeq4e`

Once the global displacement vector is solved, we extract the local displacement vector $\mathbf{d}_e$ and use it to compute the internal forces, stresses, and strains at the Gauss points using `qeq4e`.

```python
from femlabpy.elements.quads import qeq4e

# Assume some arbitrary computed displacements for the 8 DOFs
Ue = np.array([0.0, 0.0, 1e-4, -5e-5, 1.2e-4, -2e-5, 0.2e-4, 1e-5])

# Compute internal forces (qe), stresses (Se), and strains (Ee)
qe, Se, Ee = qeq4e(Xe, Ge, Ue)

print("Internal Force Vector (qe):\n", qe.ravel())
print("\nStresses at 4 Gauss Points (Sxx, Syy, Sxy):\n", Se)
```

The function `qeq4e` evaluates $\boldsymbol{\varepsilon} = \mathbf{B}\mathbf{u}_e$ and $\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$ at each of the four Gauss points. It also accumulates the internal force vector $\mathbf{q}_e = \int \mathbf{B}^T \boldsymbol{\sigma} |\mathbf{J}| t \, d\xi d\eta$.

### 5.3 Scalar Potential Flow: `keq4p`

The Q4 formulation easily adapts to scalar field problems governed by the Poisson or Laplace equations, such as steady-state heat conduction, groundwater seepage, or ideal fluid flow. In these problems, there is only one degree of freedom (potential, e.g., temperature) per node.

The "strain" is replaced by the potential gradient $\nabla u$, and the constitutive matrix is replaced by the conductivity tensor $\mathbf{k}$. The element conductivity matrix (analogous to the stiffness matrix) is $4 \times 4$.

```python
from femlabpy.elements.quads import keq4p

# Material property for potential flow: [k, b]
# k: isotropic conductivity, b: volumetric reaction term (optional)
Ge_thermal = np.array([45.0, 0.0]) # Example: thermal conductivity of steel

# Compute the 4x4 conductivity matrix
Ke_thermal = keq4p(Xe, Ge_thermal)
print("Thermal Conductivity Matrix Shape:", Ke_thermal.shape)
```

The `keq4p` function evaluates the integral $\mathbf{K}_e = \int \mathbf{B}^T \mathbf{k} \mathbf{B} |\mathbf{J}| \, d\xi d\eta$, where $\mathbf{B}$ now contains only the Cartesian derivatives of the shape functions without the expansion for $x$ and $y$ vectors.

## 6. Summary

The 2D Isoparametric Q4 element is a workhorse in computational mechanics. By mapping complex geometries to a simple parent domain, evaluating Jacobian matrices to relate local and global spaces, and employing Gauss-Legendre quadrature, it provides a robust numerical framework. The `femlabpy` library implements these steps efficiently, providing a transparent and rigorous toolset for engineers simulating complex 2D structural and thermal systems.
