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

The code uses a direct-elimination strategy with a scaled diagonal entry. That is close to the textbook "large spring" idea, but the implementation is more careful: before a constrained row is overwritten, its coupling terms are moved to the right-hand side so non-zero prescribed values are handled correctly.

The active helper is `setbc(K, p, C, dof)`. The constraint table `C` follows the same legacy layout used throughout the repository:

`[node, value]` for scalar problems.
`[node, local_dof, value]` for vector problems.

### How `setbc` works

The implementation in `src/femlabpy/boundary.py` does four things in order:

1. It computes a diagonal scale `ks` from the current stiffness matrix.
2. It converts the 1-based node/component pairs into global DOF indices.
3. It subtracts the coupling contribution of each constrained DOF from the global load vector.
4. It zeros the row and column, then places `ks` on the diagonal and `ks * value` in the load vector.

```python
    ks = 0.1 * max_abs_diagonal(K)
    if ks == 0.0:
        ks = 1.0

    # ... (DOF parsing omitted for brevity) ...

    for k in range(len(cdofs)):
        j = int(cdofs[k])
        val = cvals[k]
        # Transfer coupling forces to RHS *before* zeroing the column.
        if val != 0.0:
            if sparse:
                col_j = np.asarray(K[:, j].toarray()).ravel()
            else:
                col_j = K[:, j].copy()
            p[:, 0] -= col_j * val
            
        # Zero row and column, set diagonal spring.
        K[j, :] = 0
        K[:, j] = 0
        K[j, j] = ks
        p[j, 0] = ks * val
```

This procedure keeps the constrained value exactly visible in the modified system and avoids leaving stale coupling terms in the matrix.

---

## 2.3 General constraints

Some problems need multi-point constraints instead of single-node fixities. In that case the condition has the form:

$$ \mathbf{G} \mathbf{u} = \mathbf{Q} $$

where $\mathbf{G}$ is the constraint matrix and $\mathbf{Q}$ is the constraint right-hand side. `femlabpy` handles this with Lagrange multipliers through `solve_lag_general`.

### Building the saddle-point system

The augmented system couples the unknown displacements and the multipliers:

$$ \begin{bmatrix} \mathbf{K} & \mathbf{G}^T \\ \mathbf{G} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \mathbf{u} \\ \boldsymbol{\lambda} \end{bmatrix} = \begin{bmatrix} \mathbf{P} \\ \mathbf{Q} \end{bmatrix} $$

In `solve_lag_general`, the constraint rows are scaled to stay numerically compatible with the stiffness matrix. The function then builds the block matrix with `numpy.block` for dense systems, or `scipy.sparse.bmat` for sparse systems.

```python
    Gbar = scale * constraint_matrix
    Qbar = scale * constraint_rhs
    
    # Building the Saddle-Point Matrix using np.block
    Kbar = np.block(
        [
            [as_float_array(K), Gbar.T],
            [
                Gbar,
                np.zeros(
                    (constraint_matrix.shape[0], constraint_matrix.shape[0]),
                    dtype=float,
                ),
            ],
        ]
    )
```

The returned solution contains the physical displacement vector. If requested, the recovered multipliers are rescaled back to the original constraint units.

### Example: two-spring system

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
