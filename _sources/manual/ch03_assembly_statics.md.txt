# Chapter 3: Assembly & Constraints

This chapter details the numerical procedures used to assemble global systems of equations and enforce kinematic constraints.

## 3.1 Global Assembly Algorithm

In the finite element method, the global stiffness matrix $\mathbf{K}$ is constructed by summing the contributions of individual element stiffness matrices $\mathbf{K}_e$. Mathematically, this is expressed using Boolean connectivity matrices $\mathbf{L}_e$:

$$ \mathbf{K} = \sum_{e=1}^{nel} \mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e $$

In `femlabpy`, this assembly is achieved through the driver functions (e.g., `kq4e`, `kt3e`) which internally rely on indexing mechanisms. For example, for a 4-node quad element with 2 DOFs per node, the global DOF indices for the element are computed as:

```python
global_dofs = [2*n-2, 2*n-1 for n in element_nodes]
```

The $8 \times 8$ element matrix $\mathbf{K}_e$ is then added into the $N_{dof} \times N_{dof}$ global matrix at the intersection of `global_dofs` rows and columns. Similarly, internal element force vectors $\mathbf{q}_e$ are assembled into the global internal force vector $\mathbf{q}$.

## 3.2 Dirichlet Boundary Conditions

To solve the equilibrium equations $\mathbf{K} \mathbf{u} = \mathbf{p}$, the global stiffness matrix must be rendered non-singular by preventing rigid body motions. This is achieved by prescribing known displacements (Dirichlet boundary conditions).

`femlabpy` utilizes a direct modification approach (often conceptually similar to the penalty method) to enforce $u_i = \bar{u}_i$. The `setbc` function applies this using a massive artificial stiffness:

1. A very large stiffness value $k_{bc}$ is determined based on the maximum diagonal entry of $\mathbf{K}$:
   $$ k_{bc} = 10^6 \times \max(\text{diag}(\mathbf{K})) $$
2. For each constrained degree of freedom $i$, the corresponding row and column in $\mathbf{K}$ are zeroed out.
3. The diagonal entry $\mathbf{K}_{ii}$ is replaced with $k_{bc}$.
4. The load vector entry $\mathbf{p}_i$ is modified to $k_{bc} \times \bar{u}_i$.

$$
\mathbf{K} = \begin{bmatrix}
\ddots & 0 & \dots \\
0 & k_{bc} & 0 \\
\vdots & 0 & \ddots
\end{bmatrix},
\quad
\mathbf{p} = \begin{bmatrix}
\vdots \\
k_{bc} \bar{u}_i \\
\vdots
\end{bmatrix}
$$

When the system is solved, the equation for row $i$ yields $k_{bc} u_i = k_{bc} \bar{u}_i$, strictly enforcing the constraint.

## 3.3 General Linear Constraints (Lagrange Multipliers)

For advanced constraints where multiple DOFs are coupled linearly (e.g., rigid links or periodic boundaries), `femlabpy` employs the method of Lagrange multipliers. The constraint equation is defined as:

$$ \mathbf{G} \mathbf{u} = \mathbf{Q} $$

The potential energy functional is augmented with Lagrange multipliers $\lambda$, representing the constraint forces. This leads to an expanded saddle-point system:

$$
\begin{bmatrix}
\mathbf{K} & \mathbf{G}^T \\
\mathbf{G} & \mathbf{0}
\end{bmatrix}
\begin{Bmatrix}
\mathbf{u} \\
\lambda
\end{Bmatrix}
=
\begin{Bmatrix}
\mathbf{p} \\
\mathbf{Q}
\end{Bmatrix}
$$

The `solve_lag_general` function solves this indefinite system directly, yielding both the constrained displacements $\mathbf{u}$ and the precise constraint forces $\lambda$.
