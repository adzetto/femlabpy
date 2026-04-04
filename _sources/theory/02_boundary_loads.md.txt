# Chapter 2: Boundary Conditions and Loads

In computational mechanics, the formulation of a finite element problem ultimately leads to a system of algebraic equations, typically expressed in the form of a global system:

$$ \mathbf{K} \mathbf{u} = \mathbf{P} $$

where $\mathbf{K}$ is the global stiffness matrix, $\mathbf{u}$ is the vector of unknown nodal displacements, and $\mathbf{P}$ is the global load vector. However, before this system can be solved, the equations are singular; they represent a body floating freely in space with rigid body modes. To render the system nonsingular and to reflect the physical reality of the problem, we must impose constraints (Boundary Conditions) and apply external forces (Loads).

This chapter discusses the mathematical theory and the computational implementation of boundary conditions and loads within the `femlabpy` framework, specifically focusing on `src/femlabpy/boundary.py` and `src/femlabpy/loads.py`.

---

## 2.1 Applying External Loads

The global load vector $\mathbf{P}$ accumulates all external forces applied to the structure. These forces can originate from nodal point loads, surface tractions (distributed loads), or body forces (like gravity).

In standard finite element assembly, the contributions from elements are transformed into equivalent nodal loads and added into the global vector. When dealing with point loads explicitly applied to specific degrees of freedom (DOFs), we directly manipulate the entries of $\mathbf{P}$.

### 2.1.1 Setting and Adding Loads

The `femlabpy` library provides two primary functions in `loads.py` for manipulating the global load vector: `setload` and `addload`.

*   **`setload(P, dofs, values)`**: This function explicitly overrides the current value at the specified DOFs with new values. This is mathematically equivalent to setting $P_i = v_i$. It is typically used for initial load assignments or when a specific DOF must have a precise load value, ignoring previous accumulations.
*   **`addload(P, dofs, values)`**: In linear superposition, loads are cumulative. The `addload` function adds new load values to the existing values in the load vector. Mathematically, $P_i \leftarrow P_i + v_i$. This is the standard method for assembling structural loads from multiple overlapping sources.

```python
# Example: Applying loads in femlabpy
import numpy as np
from femlabpy.loads import setload, addload

# Initialize a global load vector for a 6-DOF system
P = np.zeros(6)

# Set a specific load at DOF 2 to 50.0
setload(P, [2], [50.0])

# Add multiple loads: 10.0 at DOF 2, and -25.0 at DOF 4
addload(P, [2, 4], [10.0, -25.0])

print("Global Load Vector:", P)
# Expected Output: [0.  0.  60. 0. -25. 0.]
```

---

## 2.2 Dirichlet Boundary Conditions: The Penalty Method

Dirichlet boundary conditions, also known as essential boundary conditions, specify the known values of the primary field variable (e.g., displacements) at certain boundaries. Mathematically, for a set of constrained DOFs $\mathcal{C}$, we require:

$$ u_i = \bar{u}_i \quad \forall i \in \mathcal{C} $$

Enforcing this directly in the system $\mathbf{K} \mathbf{u} = \mathbf{P}$ requires modifying the equations. While condensation (row/column elimination) is a mathematically exact approach, it alters the size of the matrix, making bookkeeping complex. A highly efficient alternative heavily utilized in structural analysis is the **Penalty Method** (or "Large Spring" approach).

### 2.2.1 Mathematics of the Large Spring Method

The core idea is to attach an extremely stiff spring to the constrained DOF and apply an enormous force such that the spring extends exactly by the desired displacement $\bar{u}_i$. 

Let the maximum diagonal element of the stiffness matrix be $K_{max} = \max_{j} (K_{jj})$. We define a penalty parameter $\alpha$ that is significantly larger than the stiffness of the system, typically:

$$ \alpha = 10^6 \times K_{max} $$

To enforce $u_i = \bar{u}_i$, we modify the system as follows:
1.  **Modify Stiffness:** Add $\alpha$ to the diagonal entry of $\mathbf{K}$:
    $$ K_{ii} \leftarrow K_{ii} + \alpha \approx \alpha $$
2.  **Modify Load:** Set the corresponding entry in the load vector to:
    $$ P_i \leftarrow \alpha \times \bar{u}_i $$

Since $\alpha$ is overwhelmingly large compared to other terms in the $i$-th equation, the equation effectively reduces to:

$$ \alpha u_i = \alpha \bar{u}_i \implies u_i \approx \bar{u}_i $$

This approach preserves the symmetry, size, and band structure of $\mathbf{K}$, making it ideal for sparse matrix solvers.

### 2.2.2 Implementation in `femlabpy`

The `setbc` function in `src/femlabpy/boundary.py` implements this Large Spring method perfectly.

```python
# Example: Enforcing Dirichlet BCs using the Penalty Method
import numpy as np
from femlabpy.boundary import setbc

# Define a hypothetical 3x3 stiffness matrix and load vector
K = np.array([[1000.0, -500.0, 0.0],
              [-500.0, 1000.0, -500.0],
              [0.0, -500.0, 500.0]])
P = np.array([0.0, 100.0, 50.0])

# We want to constrain DOF 0 to 0.0 (fixed support)
dofs_to_constrain = [0]
prescribed_values = [0.0]

# Apply the boundary conditions in place
setbc(K, P, dofs_to_constrain, prescribed_values)

print("Modified K_{0,0}:", K[0,0]) 
# Will be 1000.0 + 1e6 * 1000.0 = 1000001000.0
```

---

## 2.3 General Constraints: Lagrange Multipliers

While Dirichlet conditions fix individual DOFs, many engineering problems involve multi-point constraints (MPCs), where a linear combination of DOFs must satisfy a condition:

$$ c_1 u_1 + c_2 u_2 + \dots + c_n u_n = Q $$

In matrix form, a set of $m$ general constraints on $n$ DOFs can be expressed as:

$$ \mathbf{G} \mathbf{u} = \mathbf{Q} $$

where $\mathbf{G}$ is an $m \times n$ constraint matrix and $\mathbf{Q}$ is an $m \times 1$ vector of constraint values. The penalty method struggles with complex MPCs because it couples terms indiscriminately, potentially destroying matrix conditioning. Instead, we use the method of **Lagrange Multipliers**.

### 2.3.1 Augmented System Formulation

We introduce a vector of Lagrange multipliers, $\boldsymbol{\lambda}$ (of size $m \times 1$), which physically represent the constraint forces required to maintain the specified relations. The total potential energy functional is augmented to include the constraint:

$$ \Pi^*(\mathbf{u}, \boldsymbol{\lambda}) = \frac{1}{2} \mathbf{u}^T \mathbf{K} \mathbf{u} - \mathbf{u}^T \mathbf{P} + \boldsymbol{\lambda}^T (\mathbf{G} \mathbf{u} - \mathbf{Q}) $$

Taking the variation with respect to both $\mathbf{u}$ and $\boldsymbol{\lambda}$ yields the saddle-point equations:

$$ \begin{bmatrix} \mathbf{K} & \mathbf{G}^T \\ \mathbf{G} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \mathbf{u} \\ \boldsymbol{\lambda} \end{bmatrix} = \begin{bmatrix} \mathbf{P} \\ \mathbf{Q} \end{bmatrix} $$

Solving this augmented system gives both the exact displacements $\mathbf{u}$ satisfying the constraints and the reaction forces $\boldsymbol{\lambda}$ simultaneously. Note that the augmented matrix loses positive-definiteness (it has zero blocks on the diagonal), requiring robust solvers (like LDLT or LU decomposition).

### 2.3.2 Using `solve_lag_general`

The `femlabpy` framework handles this via the `solve_lag_general` function.

```python
# Example: Using Lagrange Multipliers for an MPC
import numpy as np
from femlabpy.boundary import solve_lag_general

# Base system
K = np.array([[200., -100.], [-100., 100.]])
P = np.array([0., 50.])

# Constraint: u_0 - u_1 = 0.5 (Nodes move with a relative offset of 0.5)
G = np.array([[1.0, -1.0]])
Q = np.array([0.5])

# Solve the augmented system
u, lambda_vec = solve_lag_general(K, P, G, Q)

print("Displacements:", u)
print("Constraint Force (Lagrange Multiplier):", lambda_vec)
```

---

## 2.4 Extraction of Reaction Forces

After solving the global system $\mathbf{K} \mathbf{u} = \mathbf{P}$ subject to boundary conditions, a critical post-processing step is evaluating the reaction forces at the supports. These represent the forces the ground (or constraints) exert back onto the structure.

### 2.4.1 Equilibrium and Reactions

Before constraints are applied, equilibrium states that internal forces $\mathbf{F}_{int} = \mathbf{K}_{original} \mathbf{u}$ must balance external forces $\mathbf{P}_{applied}$ plus reaction forces $\mathbf{R}$:

$$ \mathbf{K}_{original} \mathbf{u} = \mathbf{P}_{applied} + \mathbf{R} $$

Rearranging this gives the fundamental equation for reaction extraction:

$$ \mathbf{R} = \mathbf{K}_{original} \mathbf{u} - \mathbf{P}_{applied} $$

**Crucial Note:** In this formula, $\mathbf{K}_{original}$ and $\mathbf{P}_{applied}$ MUST be the *unmodified* matrices before the penalty method (`setbc`) altered them. 

For a DOF $i$ that is unconstrained, the equation evaluates to $R_i \approx 0$ (modulo numerical roundoff), verifying equilibrium. For a constrained DOF, $R_i$ will equal the reaction force required to maintain $\bar{u}_i$.

### 2.4.2 Implementation: The `reaction` Function

In `femlabpy`, the `reaction` function automates this computation, returning the reaction forces specifically at the requested DOFs.

```python
# Example: Extracting Reaction Forces
import numpy as np
from femlabpy.loads import reaction

# Assuming K_orig and P_orig were saved BEFORE applying boundary conditions
# and 'u' is the solved displacement vector.
K_orig = np.array([[1000., -1000.], [-1000., 1000.]])
P_orig = np.array([0., 500.])
u = np.array([0.0, 0.5]) # Known from solving (Node 0 fixed, Node 1 pulled)

# Extract reactions at the fixed DOF (node 0)
fixed_dofs = [0]
R = reaction(K_orig, u, P_orig, fixed_dofs)

print(f"Reaction force at DOF 0: {R[0]}")
# mathematically: R_0 = (1000*0 + -1000*0.5) - 0 = -500.
```

## Summary

Proper handling of boundary conditions and external loads dictates the stability and correctness of finite element solutions. `femlabpy` combines the computationally efficient Penalty Method for standard Dirichlet conditions with the rigorous Lagrange Multiplier formulation for complex multi-point constraints. Utilizing `setload`, `addload`, `setbc`, `solve_lag_general`, and `reaction` forms the complete cycle of setting up and recovering boundary data in computational mechanics models.