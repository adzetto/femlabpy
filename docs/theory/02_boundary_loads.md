# Boundary conditions and loads

In `femlabpy`, a finite element model eventually becomes a global linear system with known and unknown degrees of freedom. The solver needs two things before that system can be used directly: the external load vector and the boundary-condition handling that removes or constrains rigid body motion.

$$ \mathbf{K} \mathbf{u} = \mathbf{P} $$

where $\mathbf{K}$ is the global stiffness matrix, $\mathbf{u}$ is the vector of nodal unknowns, and $\mathbf{P}$ is the global load vector. Before the system is solved, the unconstrained equations are usually singular, so the code must first apply loads and then enforce constraints in a way that matches the stored DOF layout.

This chapter follows the actual implementation in `loads.py` and `boundary.py`. The important point is that the public helpers do not hide the array shape: they expect the same column conventions used everywhere else in the codebase, and the theory here explains why.

---

## 2.1 Applying external loads

The global load vector $\mathbf{P}$ accumulates all external forces applied to the structure. In this repository the most common case is a nodal load table where each row begins with a 1-based node index and the remaining columns are the force components in DOF order.

The helpers in `src/femlabpy/loads.py` do two different things:

`setload(p, P)` replaces the current values at the listed DOFs.
`addload(p, P)` accumulates additional values at the same DOFs.

That distinction matters in assembly workflows. Use `setload` when the input table is the full load definition for the problem. Use `addload` when several tables or source terms need to be combined into one vector.

### Load table layout

```python
import numpy as np
from femlabpy import init
from femlabpy.loads import setload, addload

_, p, _ = init(nn=3, dof=2)

# [node, Fx, Fy]
P = np.array([
    [1, 0.0, -100.0],
    [3, 50.0, 0.0],
])

p = setload(p, P)
p = addload(p, np.array([[3, 0.0, -25.0]]))
```

---

## 2.2 Dirichlet boundary conditions

Dirichlet boundary conditions specify known values of the primary field variable at certain nodes or components:

$$ u_i = \bar{u}_i \quad \forall i \in \mathcal{C} $$

### Mathematical formulation

Consider the partitioned system:

$$
\begin{bmatrix} \mathbf{K}_{ff} & \mathbf{K}_{fc} \\ \mathbf{K}_{cf} & \mathbf{K}_{cc} \end{bmatrix} \begin{bmatrix} \mathbf{u}_f \\ \mathbf{u}_c \end{bmatrix} = \begin{bmatrix} \mathbf{P}_f \\ \mathbf{P}_c \end{bmatrix}
$$

where subscript $f$ denotes free DOFs and $c$ denotes constrained DOFs with $\mathbf{u}_c = \bar{\mathbf{u}}_c$.

**Condensed system:** The free DOFs satisfy:

$$
\mathbf{K}_{ff} \mathbf{u}_f = \mathbf{P}_f - \mathbf{K}_{fc} \bar{\mathbf{u}}_c
$$

The term $\mathbf{K}_{fc} \bar{\mathbf{u}}_c$ transfers the effect of prescribed displacements to the right-hand side.

**Reaction forces:** After solving for $\mathbf{u}_f$:

$$
\mathbf{R}_c = \mathbf{K}_{cf} \mathbf{u}_f + \mathbf{K}_{cc} \bar{\mathbf{u}}_c - \mathbf{P}_c
$$

### Penalty method alternative

Instead of elimination, constrained DOFs can be enforced with a large spring:

$$
K_{ii} \leftarrow K_{ii} + k_s, \quad P_i \leftarrow P_i + k_s \bar{u}_i
$$

where $k_s \gg \max(K_{ii})$. This approximately enforces $u_i \approx \bar{u}_i$.

### Implementation in `setbc`

The code uses a hybrid approach that maintains matrix symmetry:

1. Compute scale: $k_s = 0.1 \cdot \max_i |K_{ii}|$
2. Transfer coupling to RHS: $\mathbf{P} \leftarrow \mathbf{P} - \mathbf{K}_{:,j} \bar{u}_j$
3. Zero row and column: $K_{j,:} = 0$, $K_{:,j} = 0$
4. Set diagonal and load: $K_{jj} = k_s$, $P_j = k_s \bar{u}_j$

This preserves symmetry and gives exact constraint enforcement.

---

## 2.3 Multi-point constraints (Lagrange multipliers)

Multi-point constraints couple multiple DOFs:

$$ \mathbf{G} \mathbf{u} = \mathbf{Q} $$

where $\mathbf{G}$ is the $m \times n$ constraint matrix and $\mathbf{Q}$ is the $m \times 1$ constraint RHS.

### Variational formulation

The constrained minimization problem:

$$
\min_{\mathbf{u}} \frac{1}{2}\mathbf{u}^T\mathbf{K}\mathbf{u} - \mathbf{u}^T\mathbf{P} \quad \text{subject to} \quad \mathbf{G}\mathbf{u} = \mathbf{Q}
$$

The Lagrangian is:

$$
\mathcal{L}(\mathbf{u}, \boldsymbol{\lambda}) = \frac{1}{2}\mathbf{u}^T\mathbf{K}\mathbf{u} - \mathbf{u}^T\mathbf{P} + \boldsymbol{\lambda}^T(\mathbf{G}\mathbf{u} - \mathbf{Q})
$$

**Stationarity conditions:**

$$
\frac{\partial \mathcal{L}}{\partial \mathbf{u}} = \mathbf{K}\mathbf{u} - \mathbf{P} + \mathbf{G}^T\boldsymbol{\lambda} = \mathbf{0}
$$

$$
\frac{\partial \mathcal{L}}{\partial \boldsymbol{\lambda}} = \mathbf{G}\mathbf{u} - \mathbf{Q} = \mathbf{0}
$$

### Saddle-point system

The KKT conditions form the saddle-point system:

$$ \begin{bmatrix} \mathbf{K} & \mathbf{G}^T \\ \mathbf{G} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \mathbf{u} \\ \boldsymbol{\lambda} \end{bmatrix} = \begin{bmatrix} \mathbf{P} \\ \mathbf{Q} \end{bmatrix} $$

**Physical interpretation of $\boldsymbol{\lambda}$:**

The Lagrange multipliers represent the constraint forces. From the first equation:

$$
\mathbf{K}\mathbf{u} = \mathbf{P} - \mathbf{G}^T\boldsymbol{\lambda}
$$

The term $-\mathbf{G}^T\boldsymbol{\lambda}$ is the reaction force required to enforce the constraint.

### Numerical conditioning

The saddle-point matrix is indefinite. For numerical stability, `solve_lag_general` scales the constraint equations:

$$
\bar{\mathbf{G}} = s \cdot \mathbf{G}, \quad \bar{\mathbf{Q}} = s \cdot \mathbf{Q}
$$

where $s \sim 0.01 \cdot \max_i |K_{ii}|$ balances the constraint rows with the stiffness entries.

After solving, the multipliers are rescaled: $\boldsymbol{\lambda} = s \cdot \bar{\boldsymbol{\lambda}}$

### Example: rigid link constraint

This example shows the same structure as the implementation: one equilibrium solve, one constraint equation, and one multiplier that reports the transmitted constraint force.

```python
import numpy as np

def solve_lag_general(K, p, G, Q):
    scale = 1.0e-2 * np.max(np.abs(np.diag(K)))
    Gbar = scale * G
    Qbar = scale * Q
    
    # Build the saddle-point matrix
    Kbar = np.block([
        [K, Gbar.T],
        [Gbar, np.zeros((G.shape[0], G.shape[0]))]
    ])
    
    # Build the augmented right-hand side
    pbar = np.vstack([p, Qbar])
    
    # Solve the augmented system
    augmented = np.linalg.solve(Kbar, pbar)
    u = augmented[:K.shape[0]]
    lagrange = augmented[K.shape[0]:] * scale
    return u, lagrange

def main():
    K = np.array([
        [1500.0, -500.0],
        [-500.0,  500.0]
    ])
    
    p = np.array([[0.0], 
                  [100.0]])

    G = np.array([[1.0, -1.0]])
    Q = np.array([[0.0]])

    u, lagrange = solve_lag_general(K, p, G, Q)
    
    print("Displacements:")
    print(f"u_1 = {u[0,0]:.5f}")
    print(f"u_2 = {u[1,0]:.5f}")
    print("\nConstraint Force (Lagrange Multiplier):")
    print(f"lambda = {lagrange[0,0]:.5f}")

if __name__ == "__main__":
    main()
```

The key outcome is that the constraint is enforced without guessing a penalty stiffness. That makes the method useful whenever the exact multiplier value matters, such as in periodic constraints or multi-point tie conditions.

---

## Summary

`femlabpy` uses direct elimination for standard fixed DOFs and Lagrange multipliers for general linear constraints. The two paths share the same global indexing rules, so the boundary-condition code stays consistent with the rest of the assembly pipeline.
