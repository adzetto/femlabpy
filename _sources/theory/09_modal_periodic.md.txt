# Chapter 9: Eigenvalue Problems and Periodic Boundaries

In computational mechanics, evaluating the dynamic characteristics and effective macroscopic properties of materials forms a cornerstone of modern engineering analysis. This chapter delves into two advanced topics fundamental to finite element analysis (FEA): Free Vibration Modal Analysis and the enforcement of Periodic Boundary Conditions (PBCs) for Computational Homogenization.

By examining the generalized eigenvalue problem, we will characterize the natural frequencies and mode shapes of a structure. Subsequently, we explore Representative Volume Elements (RVEs) and how microstructural properties are scaled up using periodic boundaries and Lagrange multipliers to extract the effective macroscopic matrix $\mathbf{C}_{eff}$.

---

## 9.1 Free Vibration Modal Analysis

The dynamic response of an undamped structural system without external forcing is governed by the homogeneous equation of motion:

$$
\mathbf{M} \ddot{\mathbf{u}}(t) + \mathbf{K} \mathbf{u}(t) = \mathbf{0}
$$

where $\mathbf{M}$ is the global mass matrix, $\mathbf{K}$ is the global stiffness matrix, and $\mathbf{u}(t)$ is the nodal displacement vector as a function of time. Assuming a harmonic solution of the form $\mathbf{u}(t) = \phi e^{i\omega t}$, taking the second time derivative yields $\ddot{\mathbf{u}}(t) = -\omega^2 \phi e^{i\omega t}$. Substituting this back into the equation of motion provides the classic generalized eigenvalue problem:

$$
\left( \mathbf{K} - \omega^2 \mathbf{M} \right) \phi = \mathbf{0}
$$

Rearranging this expression yields the canonical form for free vibration modal analysis:

$$
\mathbf{K}\phi = \omega^2 \mathbf{M}\phi
$$

In this formulation:
- $\omega_i^2$ represents the $i$-th eigenvalue, where $\omega_i$ is the natural angular frequency of the $i$-th mode in radians per second.
- $\phi_i$ is the $i$-th eigenvector, which physically represents the mode shape or spatial deformation pattern of the structure oscillating at $\omega_i$.

Because $\mathbf{K}$ and $\mathbf{M}$ are typically real, symmetric, and positive-definite (assuming proper rigid body constraints are applied), all eigenvalues $\omega_i^2$ are real and non-negative, and the corresponding eigenvectors are real-valued.

### 9.1.1 Mass-Normalization of Mode Shapes

Eigenvectors only represent a shape, not an absolute magnitude. For mathematical convenience and numerical stability, it is standard practice to normalize these mode shapes with respect to the mass matrix. The mass-orthogonality property of eigenvectors states that modes associated with different natural frequencies are orthogonal with respect to the mass and stiffness matrices. 

When mode shapes are mass-normalized, the generalized modal mass matrix becomes the identity matrix $\mathbf{I}$:

$$
\mathbf{\Phi}^T \mathbf{M} \mathbf{\Phi} = \mathbf{I}
$$

where $\mathbf{\Phi} = \left[ \phi_1, \phi_2, \dots, \phi_n \right]$ is the modal matrix containing the mode shapes as columns. For a single mode $i$, this reduces to:

$$
\phi_i^T \mathbf{M} \phi_i = 1
$$

Consequently, the stiffness matrix is diagonalized to contain the squared natural frequencies:

$$
\mathbf{\Phi}^T \mathbf{K} \mathbf{\Phi} = \mathbf{\Omega}^2 = \text{diag}(\omega_1^2, \omega_2^2, \dots, \omega_n^2)
$$

### 9.1.2 Modal Participation Factors (MPF)

To understand how readily a particular mode is excited by a uniform base acceleration (e.g., in earthquake engineering), we calculate the Modal Participation Factor (MPF), denoted as $\Gamma_i$. The MPF in a specific spatial direction $j$ (where $j \in \{x, y, z\}$) is defined as:

$$
\Gamma_{ij} = \frac{\phi_i^T \mathbf{M} \mathbf{r}_j}{\phi_i^T \mathbf{M} \phi_i}
$$

Since we have mass-normalized the modes ($\phi_i^T \mathbf{M} \phi_i = 1$), the equation simplifies to:

$$
\Gamma_{ij} = \phi_i^T \mathbf{M} \mathbf{r}_j
$$

Here, $\mathbf{r}_j$ is the influence vector representing the displacements of the structural degrees of freedom resulting from a unit rigid body displacement in the $j$-th direction. 

### 9.1.3 Effective Modal Mass

The Effective Modal Mass (EMM) measures the amount of system mass participating in a given mode for a specific excitation direction. It is a critical metric for determining how many modes are required to capture the bulk dynamic response of the system (typically, standards require capturing at least 90% of the total mass).

The effective modal mass for the $i$-th mode in the $j$-th direction is defined as:

$$
m_{eff, ij} = \frac{(\phi_i^T \mathbf{M} \mathbf{r}_j)^2}{\phi_i^T \mathbf{M} \phi_i}
$$

Using mass-normalized eigenvectors, this simply becomes the square of the participation factor:

$$
m_{eff, ij} = \Gamma_{ij}^2
$$

The sum of the effective modal masses for all modes equals the total mass of the structure.

```python
# Snippet demonstrating modal analysis using femlabpy
import numpy as np
from femlabpy.modal import solve_eigenproblem, compute_mpf_emm

# Assuming global stiffness K and mass M matrices are built
# Solve for first 10 modes
eigenvalues, eigenvectors = solve_eigenproblem(K, M, num_modes=10)

# Frequencies in Hz
frequencies = np.sqrt(eigenvalues) / (2 * np.pi)

# Influence vector for X-direction
r_x = np.zeros(K.shape[0])
r_x[0::2] = 1.0  # Assuming 2D, alternating X/Y DOFs

# Compute Modal Participation Factors and Effective Modal Mass
MPF_x, EMM_x = compute_mpf_emm(eigenvectors, M, r_x)

print(f"Mode 1 Frequency: {frequencies[0]:.2f} Hz")
print(f"Mode 1 Effective Mass (X): {EMM_x[0]:.2f}")
```

---

## 9.2 Periodic Boundary Conditions (PBC)

In multiscale computational mechanics, the macroscopic behavior of heterogeneous materials (like composites, foams, or lattice structures) is derived by analyzing a microscopic Representative Volume Element (RVE). An RVE must be statistically representative of the bulk material.

To ensure that the RVE behaves as though it is embedded within an infinite periodic medium, Periodic Boundary Conditions (PBCs) are applied to its boundaries. This avoids artificial stiffening or softening at the edges, which would occur with standard Dirichlet or Neumann boundaries.

### 9.2.1 The Periodic Constraint Equation

Consider a 2D rectangular RVE with dimensions $\Delta x$ and $\Delta y$. The boundaries consist of paired sets of nodes: a right face (+) and a left face (-), as well as a top face (+) and a bottom face (-).

For the RVE to remain continuous with its periodic neighbors under a given macroscopic strain field $\bar{\boldsymbol{\epsilon}}$, the displacements of corresponding nodes on opposite boundary faces must satisfy:

$$
\mathbf{u}^+ - \mathbf{u}^- = \bar{\boldsymbol{\epsilon}} \Delta \mathbf{x}
$$

For a 2D continuum, the macroscopic strain tensor has three independent components: $\bar{\epsilon}_{xx}$, $\bar{\epsilon}_{yy}$, and $\bar{\gamma}_{xy}$. The constraint equations for the $x$-displacements ($u$) and $y$-displacements ($v$) between a node on the right face $(x^+, y)$ and a corresponding node on the left face $(x^-, y)$ become:

$$
u^+ - u^- = \bar{\epsilon}_{xx} \Delta x + \frac{1}{2} \bar{\gamma}_{xy} \Delta y
$$
$$
v^+ - v^- = \frac{1}{2} \bar{\gamma}_{xy} \Delta x + \bar{\epsilon}_{yy} \Delta y
$$

### 9.2.2 Enforcement via Lagrange Multipliers

In the finite element framework, these multipoint constraints are often enforced using Lagrange multipliers. For each pair of periodic nodes, an additional algebraic constraint equation is appended to the global system. 

Let $\mathbf{C}_{pbc}$ be the boolean constraint matrix and $\mathbf{q}$ be the vector of prescribed relative displacements (derived from $\bar{\boldsymbol{\epsilon}} \Delta \mathbf{x}$). The extended system to solve becomes a saddle-point problem:

$$
\begin{bmatrix}
\mathbf{K} & \mathbf{C}_{pbc}^T \\
\mathbf{C}_{pbc} & \mathbf{0}
\end{bmatrix}
\begin{Bmatrix}
\mathbf{u} \\
\boldsymbol{\lambda}
\end{Bmatrix}
=
\begin{Bmatrix}
\mathbf{f}_{ext} \\
\mathbf{q}
\end{Bmatrix}
$$

Here, $\boldsymbol{\lambda}$ represents the Lagrange multipliers, which physically correspond to the unknown reaction forces maintaining periodicity. For homogenization, $\mathbf{f}_{ext} = \mathbf{0}$ for internal nodes, and the deformation is entirely driven by $\mathbf{q}$.

---

## 9.3 Computational Homogenization

The primary objective of analyzing an RVE under PBCs is to extract the effective macroscopic properties—specifically, the homogenized stiffness matrix (or compliance matrix) $\mathbf{C}_{eff}$. In 2D plane elasticity, $\mathbf{C}_{eff}$ is a $3 \times 3$ matrix relating the macroscopic stresses $\bar{\boldsymbol{\sigma}}$ to macroscopic strains $\bar{\boldsymbol{\epsilon}}$:

$$
\bar{\boldsymbol{\sigma}} = \mathbf{C}_{eff} \bar{\boldsymbol{\epsilon}}
$$

### 9.3.1 Driving Macro-Strain States

To fully populate the $\mathbf{C}_{eff}$ matrix, the `homogenize` routine systematically drives the RVE through three independent, unit macroscopic strain states. The three tests are:

1.  **Test 1: Pure X-Tension:**
    $$ \bar{\boldsymbol{\epsilon}}^{(1)} = [1, 0, 0]^T \implies \bar{\epsilon}_{xx} = 1, \bar{\epsilon}_{yy} = 0, \bar{\gamma}_{xy} = 0 $$
2.  **Test 2: Pure Y-Tension:**
    $$ \bar{\boldsymbol{\epsilon}}^{(2)} = [0, 1, 0]^T \implies \bar{\epsilon}_{xx} = 0, \bar{\epsilon}_{yy} = 1, \bar{\gamma}_{xy} = 0 $$
3.  **Test 3: Pure In-Plane Shear:**
    $$ \bar{\boldsymbol{\epsilon}}^{(3)} = [0, 0, 1]^T \implies \bar{\epsilon}_{xx} = 0, \bar{\epsilon}_{yy} = 0, \bar{\gamma}_{xy} = 1 $$

### 9.3.2 Extracting the Effective Stiffness

For each load case $k \in \{1, 2, 3\}$, the system of equations is solved for the nodal displacements $\mathbf{u}^{(k)}$ and the Lagrange multipliers $\boldsymbol{\lambda}^{(k)}$.

The macroscopic stress $\bar{\boldsymbol{\sigma}}^{(k)}$ corresponding to the imposed strain $\bar{\boldsymbol{\epsilon}}^{(k)}$ is calculated by averaging the stress over the volume $V$ of the RVE. Due to Hill's energy condition, this volume average can be elegantly computed using the boundary traction forces (the Lagrange multipliers) and their moment arms, or by evaluating the total strain energy:

$$
\bar{\sigma}_{ij} = \frac{1}{V} \int_V \sigma_{ij} dV = \frac{1}{V} \sum_{n \in \partial V} f_i^{(n)} x_j^{(n)}
$$

Since the imposed strain vectors were chosen as unit basis vectors, the resulting macroscopic stress vector $\bar{\boldsymbol{\sigma}}^{(k)} = [\bar{\sigma}_{xx}^{(k)}, \bar{\sigma}_{yy}^{(k)}, \bar{\sigma}_{xy}^{(k)}]^T$ directly constitutes the $k$-th column of the homogenized matrix $\mathbf{C}_{eff}$.

$$
\mathbf{C}_{eff} = \left[ \bar{\boldsymbol{\sigma}}^{(1)} \quad \bar{\boldsymbol{\sigma}}^{(2)} \quad \bar{\boldsymbol{\sigma}}^{(3)} \right]
$$

```python
# Snippet demonstrating RVE Homogenization
from femlabpy.periodic import apply_pbc, homogenize

# RVE Geometry bounds
bounds = {'x_min': 0.0, 'x_max': 1.0, 'y_min': 0.0, 'y_max': 1.0}

# Generate constraint matrix C_pbc for the RVE
C_pbc = apply_pbc(mesh, bounds)

# Compute the effective 3x3 compliance/stiffness matrix C_eff
# Internally applies 3 macro-strain states using Lagrange multipliers
C_eff = homogenize(K_global, C_pbc, bounds)

print("Effective Homogenized Matrix C_eff:")
print(np.array_str(C_eff, precision=3))
```

### 9.3.3 Summary of the Homogenization Routine

1.  **Identify Boundary Node Pairs:** Locate nodes on opposite faces of the RVE.
2.  **Formulate Constraints:** Construct the Boolean matrix $\mathbf{C}_{pbc}$.
3.  **Loop Over Strain States:** For $k=1, 2, 3$, define unit macro-strains $\bar{\boldsymbol{\epsilon}}^{(k)}$ and compute the constraint vectors $\mathbf{q}^{(k)}$.
4.  **Solve System:** Assemble and solve the saddle-point block matrix for displacements $\mathbf{u}^{(k)}$ and multipliers $\boldsymbol{\lambda}^{(k)}$.
5.  **Compute Macro-Stress:** Integrate boundary forces or element stresses to find $\bar{\boldsymbol{\sigma}}^{(k)}$.
6.  **Assemble $\mathbf{C}_{eff}$:** Combine the column vectors into the final $3 \times 3$ effective matrix.

The resulting $\mathbf{C}_{eff}$ seamlessly replaces detailed microstructural models in macro-scale FE simulations, retaining the intricate geometric and material influences of the underlying architecture.
