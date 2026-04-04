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

# Chapter 3: 1D Bar Elements (Linear & Nonlinear)

This chapter introduces the formulation of 1D bar elements (trusses), which are fundamental building blocks in computational mechanics. We will start with the basic linear assumption and extend the theory to encompass large deformations using the Green-Lagrange strain. 

By the end of this chapter, you will understand the underlying physics and mathematics, and see exactly how the tangent stiffness matrices and internal force vectors are implemented in Python.

## 3.1 Linear 1D Bar Formulation

The simplest finite element is the 2-node linear bar element. A bar can only transmit axial loads—meaning it has no bending or shear stiffness.

### 3.1.1 Kinematics and Strain
Consider a 1D bar of initial length $L$, cross-sectional area $A$, and Young's modulus $E$. Let the bar be aligned along the local $x$-axis. The bar is defined by two nodes, $1$ and $2$, with displacements $u_1$ and $u_2$.

The axial strain $\varepsilon$ under the assumption of small deformations is simply the change in length divided by the original length:

$$ \varepsilon = \frac{\Delta L}{L} = \frac{u_2 - u_1}{L} $$

Using the standard finite element approach, we define the element displacement vector $\mathbf{d}^e = [u_1, u_2]^T$. The strain-displacement relation can be written in matrix form:

$$ \varepsilon = \mathbf{B} \mathbf{d}^e = \begin{bmatrix} -\frac{1}{L} & \frac{1}{L} \end{bmatrix} \begin{bmatrix} u_1 \\ u_2 \end{bmatrix} $$

where $\mathbf{B}$ is the strain-displacement matrix.

### 3.1.2 The Linear Stiffness Matrix
The principle of virtual work or minimizing the total potential energy yields the element stiffness matrix $\mathbf{K}^e$:

$$ \mathbf{K}^e = \int_V \mathbf{B}^T E \mathbf{B} \, dV $$

Since the area $A$ and modulus $E$ are constant, the volume integral becomes $A \cdot \int_0^L dx$:

$$ \mathbf{K}^e = E A L \left( \begin{bmatrix} -\frac{1}{L} \\ \frac{1}{L} \end{bmatrix} \begin{bmatrix} -\frac{1}{L} & \frac{1}{L} \end{bmatrix} \right) = \frac{EA}{L} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix} $$

This is the canonical local linear stiffness matrix for a 2-node bar element. 

## 3.2 Geometric Nonlinear Formulation

When displacements and rotations become large, the linear assumption ($\varepsilon = \Delta L / L$) is no longer valid. Even if the material remains linear-elastic, large rigid body rotations introduce significant geometric nonlinearity.

To capture large deformations, we must use an objective measure of strain, such as the Green-Lagrange strain.

### 3.2.1 Kinematics and Green-Lagrange Strain
Let the initial node coordinates be $\mathbf{X}_1$ and $\mathbf{X}_2$. The initial length vector is $\mathbf{a}_0 = \mathbf{X}_2 - \mathbf{X}_1$, and the initial length is $l_0 = \|\mathbf{a}_0\|$.

During deformation, the nodes move to new current coordinates $\mathbf{x}_1 = \mathbf{X}_1 + \mathbf{u}_1$ and $\mathbf{x}_2 = \mathbf{X}_2 + \mathbf{u}_2$. The current length vector is $\mathbf{a}_1 = \mathbf{x}_2 - \mathbf{x}_1$, and the current length is $l_1 = \|\mathbf{a}_1\|$.

The 1D Green-Lagrange strain $E_{GL}$ is defined by the squared change in lengths:

$$ \varepsilon = E_{GL} = \frac{l_1^2 - l_0^2}{2 l_0^2} $$

This strain measure is invariant under rigid body translations and rotations because it strictly depends on the squared lengths. 

### 3.2.2 Internal Force Response ($\mathbf{q}^e$)
The internal force vector is derived from the variation of the internal strain energy $U = \int_V \frac{1}{2} E \varepsilon^2 dV$. Using the chain rule on the strain variation $\delta \varepsilon$:

$$ \delta \varepsilon = \frac{1}{l_0^2} \mathbf{a}_1 \cdot \delta \mathbf{a}_1 = \frac{1}{l_0^2} \mathbf{a}_1^T (\delta \mathbf{u}_2 - \delta \mathbf{u}_1) $$

Defining the normal stress as $S = E \varepsilon$ and the normal force as $N = A S = A E \varepsilon$, the internal force vector corresponding to the nodal displacements is:

$$ \mathbf{q}^e = \frac{A E \varepsilon}{l_0} \begin{bmatrix} -\mathbf{a}_1 \\ \mathbf{a}_1 \end{bmatrix} $$

#### Implementation: `qebar`
Here is the exact NumPy implementation from `femlabpy` for evaluating the internal force of a nonlinear bar.

```python
import numpy as np

def qebar(Xe0, Xe1, Ge):
    """Compute the internal-force response of a single geometrically nonlinear bar."""
    initial = np.array(Xe0, dtype=float)
    current = np.array(Xe1, dtype=float)
    props = np.array(Ge, dtype=float).reshape(-1)

    a0 = (initial[1] - initial[0]).reshape(-1, 1)
    l0 = float(np.linalg.norm(a0))
    a1 = (current[1] - current[0]).reshape(-1, 1)
    l1 = float(np.linalg.norm(a1))

    A = props[0]
    E = props[1] if props.size > 1 else 1.0
    
    # Green-Lagrange Strain
    strain = 0.5 * (l1**2 - l0**2) / l0**2
    stress = E * strain
    
    # Internal Force Vector
    qe = (A * stress / l0) * np.vstack([-a1, a1])
    return qe, float(stress), float(strain)
```

### 3.2.3 Tangent Stiffness Matrix ($\mathbf{K}^e_{tan}$)
For nonlinear solvers like the Newton-Raphson method, we require the tangent stiffness matrix, which is the derivative of the internal force vector with respect to the nodal displacements:

$$ \mathbf{K}_{tan}^e = \frac{\partial \mathbf{q}^e}{\partial \mathbf{d}^e} $$

Taking the derivative of the internal force vector yields two distinct components: the material stiffness matrix $\mathbf{K}_M$ and the geometric stiffness matrix $\mathbf{K}_G$.

$$ \mathbf{K}_{tan}^e = \mathbf{K}_M + \mathbf{K}_G $$

#### 1. Material Stiffness Matrix ($\mathbf{K}_M$)
This part comes from the derivative of the strain term. It represents the axial stiffness projected into the current configuration's spatial orientation:

$$ \mathbf{K}_M = \frac{E A}{l_0^3} \begin{bmatrix} \mathbf{a}_1 \mathbf{a}_1^T & -\mathbf{a}_1 \mathbf{a}_1^T \\ -\mathbf{a}_1 \mathbf{a}_1^T & \mathbf{a}_1 \mathbf{a}_1^T \end{bmatrix} $$

Notice that $\mathbf{a}_1 \mathbf{a}_1^T$ is an outer product creating an $N_d \times N_d$ matrix corresponding to the spatial dimension of the bar.

#### 2. Geometric Stiffness Matrix ($\mathbf{K}_G$)
This part comes from the derivative of the spatial vector $\mathbf{a}_1$ while treating the force $N$ as constant. It accounts for the stiffness modification due to existing internal tension or compression:

$$ \mathbf{K}_G = \frac{N}{l_0} \begin{bmatrix} \mathbf{I} & -\mathbf{I} \\ -\mathbf{I} & \mathbf{I} \end{bmatrix} $$

where $N = A E \varepsilon$ is the internal normal force, and $\mathbf{I}$ is the identity matrix. Under compression ($N < 0$), this term reduces the total stiffness, which can naturally lead to buckling instability when $\det(\mathbf{K}_{tan}^e) \le 0$.

#### Implementation: `kebar`
Here is the exact NumPy implementation from `femlabpy` for generating the local nonlinear tangent stiffness matrix.

```python
def kebar(Xe0, Xe1, Ge):
    """Compute the tangent stiffness matrix of a geometrically nonlinear bar element."""
    initial = np.array(Xe0, dtype=float)
    current = np.array(Xe1, dtype=float)
    props = np.array(Ge, dtype=float).reshape(-1)

    a0 = (initial[1] - initial[0]).reshape(-1, 1)
    l0 = float(np.linalg.norm(a0))
    a1 = (current[1] - current[0]).reshape(-1, 1)
    l1 = float(np.linalg.norm(a1))

    A = props[0]
    E = props[1] if props.size > 1 else 1.0
    
    # Strain and Normal Force
    strain = 0.5 * (l1**2 - l0**2) / l0**2
    normal_force = A * E * strain
    
    # Matrix Components
    identity = np.eye(a0.shape[0], dtype=float)
    
    # K_M (Material Stiffness)
    K_M = (E * A / l0**3) * np.block([
        [ a1 @ a1.T, -a1 @ a1.T],
        [-a1 @ a1.T,  a1 @ a1.T]
    ])
    
    # K_G (Geometric Stiffness)
    K_G = (normal_force / l0) * np.block([
        [ identity, -identity],
        [-identity,  identity]
    ])
    
    return K_M + K_G
```

## 3.3 The Newton-Raphson Iteration Scheme

In structural mechanics, solving a nonlinear equilibrium path $\mathbf{q}(\mathbf{d}) = \mathbf{f}_{ext}$ requires iterative methods. The Newton-Raphson scheme iteratively finds the displacement increment $\Delta \mathbf{d}$ using the tangent stiffness matrix:

1. Evaluate the internal force $\mathbf{q}^{(k)}$ using the current displacements $\mathbf{d}^{(k)}$.
2. Compute the residual force $\mathbf{r}^{(k)} = \mathbf{f}_{ext} - \mathbf{q}^{(k)}$.
3. If $\|\mathbf{r}^{(k)}\| < \text{tolerance}$, the solution has converged.
4. Otherwise, evaluate the tangent stiffness matrix $\mathbf{K}_{tan}^{(k)}$.
5. Solve for the displacement increment: $\mathbf{K}_{tan}^{(k)} \Delta \mathbf{d} = \mathbf{r}^{(k)}$.
6. Update displacements: $\mathbf{d}^{(k+1)} = \mathbf{d}^{(k)} + \Delta \mathbf{d}$.
7. Repeat the process.

This robust solver mechanism leans entirely on the accurate formulation of `kebar` and `qebar` for rapid (quadratic) convergence during finite element analysis.
