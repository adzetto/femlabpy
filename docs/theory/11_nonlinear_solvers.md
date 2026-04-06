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

# Nonlinear Solvers

This chapter derives the mathematical foundations for nonlinear finite element analysis: Newton-Raphson iteration and arc-length methods. We then connect these to the `femlabpy` implementation.

## 1. The Nonlinear Problem

### 1.1 Equilibrium equation

For a structure under load, equilibrium requires:

$$
\mathbf{R}(\mathbf{u}) = \mathbf{F}^{ext} - \mathbf{F}^{int}(\mathbf{u}) = \mathbf{0}
$$

where:
- $\mathbf{u}$ is the displacement vector
- $\mathbf{F}^{ext}$ is the external load vector
- $\mathbf{F}^{int}(\mathbf{u})$ is the internal force vector
- $\mathbf{R}$ is the residual (out-of-balance) force

Nonlinearity arises from:
- **Material nonlinearity**: $\boldsymbol{\sigma} = \boldsymbol{\sigma}(\boldsymbol{\varepsilon})$ is nonlinear (plasticity, hyperelasticity)
- **Geometric nonlinearity**: strain-displacement relation is nonlinear (large deformations)

### 1.2 Linearization

The residual at $\mathbf{u} + \Delta\mathbf{u}$ is approximated by Taylor expansion:

$$
\mathbf{R}(\mathbf{u} + \Delta\mathbf{u}) \approx \mathbf{R}(\mathbf{u}) + \frac{\partial \mathbf{R}}{\partial \mathbf{u}} \Delta\mathbf{u}
$$

The tangent stiffness matrix is:

$$
\mathbf{K}_T = -\frac{\partial \mathbf{R}}{\partial \mathbf{u}} = \frac{\partial \mathbf{F}^{int}}{\partial \mathbf{u}}
$$

For equilibrium at the new state:

$$
\mathbf{R}(\mathbf{u} + \Delta\mathbf{u}) = \mathbf{0} \implies \mathbf{K}_T \Delta\mathbf{u} = \mathbf{R}(\mathbf{u})
$$

---

## 2. Newton-Raphson Method

### 2.1 Standard Newton-Raphson

The Newton-Raphson algorithm iteratively solves the nonlinear system:

**Initialize:** $\mathbf{u}^{(0)} = \mathbf{u}_n$ (previous converged state)

**Iterate** for $k = 0, 1, 2, \ldots$:

1. Compute residual:
   $$\mathbf{R}^{(k)} = \mathbf{F}^{ext} - \mathbf{F}^{int}(\mathbf{u}^{(k)})$$

2. Check convergence:
   $$\|\mathbf{R}^{(k)}\| < \text{tol} \implies \text{converged}$$

3. Compute tangent stiffness:
   $$\mathbf{K}_T^{(k)} = \frac{\partial \mathbf{F}^{int}}{\partial \mathbf{u}}\bigg|_{\mathbf{u}^{(k)}}$$

4. Solve for correction:
   $$\mathbf{K}_T^{(k)} \Delta\mathbf{u}^{(k)} = \mathbf{R}^{(k)}$$

5. Update displacement:
   $$\mathbf{u}^{(k+1)} = \mathbf{u}^{(k)} + \Delta\mathbf{u}^{(k)}$$

### 2.2 Convergence properties

Newton-Raphson has quadratic convergence near the solution:

$$
\|\mathbf{u}^{(k+1)} - \mathbf{u}^*\| \leq C \|\mathbf{u}^{(k)} - \mathbf{u}^*\|^2
$$

This means the error squares each iteration, giving rapid convergence when close to the solution.

However, convergence requires:
- A good initial guess
- Non-singular tangent matrix
- Residual in the basin of convergence

### 2.3 Modified Newton-Raphson

To reduce computation, the tangent can be held constant:

$$
\mathbf{K}_T^{(0)} \Delta\mathbf{u}^{(k)} = \mathbf{R}^{(k)} \quad \forall k
$$

This gives linear (slower) convergence but avoids reforming and refactorizing $\mathbf{K}_T$ each iteration.

---

## 3. Load Control and Incremental Analysis

### 3.1 Load incrementation

For path-dependent problems (plasticity), the load is applied in increments:

$$
\mathbf{F}^{ext}_n = \lambda_n \hat{\mathbf{F}}
$$

where $\lambda_n$ is the load factor and $\hat{\mathbf{F}}$ is the reference load.

The solution proceeds through load steps $n = 1, 2, \ldots$:

$$
\Delta\lambda_n = \lambda_n - \lambda_{n-1}
$$

### 3.2 Predictor-corrector scheme

Each load step uses:

**Predictor:** Apply load increment, compute initial displacement estimate

$$
\Delta\mathbf{u}^{(0)} = \mathbf{K}_T^{-1} (\Delta\lambda \hat{\mathbf{F}})
$$

**Corrector:** Newton iterations to equilibrium

$$
\mathbf{K}_T \delta\mathbf{u}^{(k)} = \mathbf{R}^{(k)}, \quad \mathbf{u}^{(k+1)} = \mathbf{u}^{(k)} + \delta\mathbf{u}^{(k)}
$$

---

## 4. Arc-Length Methods

### 4.1 Motivation

Standard load control fails at limit points where $\frac{d\lambda}{d\mathbf{u}} \to 0$. At such points, the load-displacement curve has a horizontal tangent and the structure exhibits snap-through or snap-back behavior.

Arc-length methods treat the load factor $\lambda$ as an additional unknown, allowing the solver to trace the complete equilibrium path including post-buckling and softening regimes.

### 4.2 Augmented system

The equilibrium equation becomes:

$$
\mathbf{R}(\mathbf{u}, \lambda) = \lambda \hat{\mathbf{F}} - \mathbf{F}^{int}(\mathbf{u}) = \mathbf{0}
$$

With $n$ displacement DOFs and one unknown $\lambda$, we have $n$ equations and $n+1$ unknowns. An additional constraint equation is needed:

$$
g(\Delta\mathbf{u}, \Delta\lambda) = 0
$$

### 4.3 Linearization

Linearizing the equilibrium equation:

$$
\mathbf{K}_T \delta\mathbf{u} - \delta\lambda \hat{\mathbf{F}} = \mathbf{R}^{(k)}
$$

The displacement correction can be decomposed:

$$
\delta\mathbf{u} = \delta\mathbf{u}_R + \delta\lambda \, \delta\mathbf{u}_F
$$

where:
- $\delta\mathbf{u}_R = \mathbf{K}_T^{-1} \mathbf{R}^{(k)}$ (residual contribution)
- $\delta\mathbf{u}_F = \mathbf{K}_T^{-1} \hat{\mathbf{F}}$ (load contribution)

### 4.4 Spherical arc-length (Riks method)

The constraint is a hypersphere in load-displacement space:

$$
g = \Delta\mathbf{u}^T \Delta\mathbf{u} + \psi^2 \Delta\lambda^2 \|\hat{\mathbf{F}}\|^2 - \Delta l^2 = 0
$$

**Derivation of $\delta\lambda$:**

Substituting $\delta\mathbf{u} = \delta\mathbf{u}_R + \delta\lambda \, \delta\mathbf{u}_F$ into the linearized constraint:

$$
(\Delta\mathbf{u} + \delta\mathbf{u})^T (\Delta\mathbf{u} + \delta\mathbf{u}) + \psi^2 (\Delta\lambda + \delta\lambda)^2 \|\hat{\mathbf{F}}\|^2 = \Delta l^2
$$

Keeping linear terms:

$$
2\Delta\mathbf{u}^T \delta\mathbf{u}_R + 2\Delta\mathbf{u}^T \delta\mathbf{u}_F \delta\lambda + 2\psi^2 \Delta\lambda \|\hat{\mathbf{F}}\|^2 \delta\lambda = \Delta l^2 - \|\Delta\mathbf{u}\|^2 - \psi^2\Delta\lambda^2\|\hat{\mathbf{F}}\|^2
$$

Solving for $\delta\lambda$:

$$
\delta\lambda = \frac{-2\Delta\mathbf{u}^T \delta\mathbf{u}_R + \Delta l^2 - \|\Delta\mathbf{u}\|^2 - \psi^2\Delta\lambda^2\|\hat{\mathbf{F}}\|^2}{2(\Delta\mathbf{u}^T \delta\mathbf{u}_F + \psi^2 \Delta\lambda \|\hat{\mathbf{F}}\|^2)}
$$

### 4.5 Cylindrical arc-length (Crisfield method)

Uses displacement-only constraint:

$$
g = \Delta\mathbf{u}^T \Delta\mathbf{u} - \Delta l^2 = 0
$$

This yields a quadratic equation in $\delta\lambda$:

$$
a_1 \delta\lambda^2 + a_2 \delta\lambda + a_3 = 0
$$

where:
- $a_1 = \delta\mathbf{u}_F^T \delta\mathbf{u}_F$
- $a_2 = 2(\Delta\mathbf{u} + \delta\mathbf{u}_R)^T \delta\mathbf{u}_F$
- $a_3 = (\Delta\mathbf{u} + \delta\mathbf{u}_R)^T (\Delta\mathbf{u} + \delta\mathbf{u}_R) - \Delta l^2$

The root is chosen to maintain path direction (positive work criterion).

### 4.6 Orthogonal residual method (femlabpy)

The `femlabpy` implementation uses orthogonality to the predictor:

$$
\delta\mathbf{u}^T \Delta\mathbf{u}^{(0)} = 0
$$

**Derivation:**

Let $\Delta\mathbf{u}^{(0)} = \mathbf{K}_T^{-1} \hat{\mathbf{F}}$ be the reference increment. The correction satisfies:

$$
(\delta\mathbf{u}_R + \delta\lambda \, \delta\mathbf{u}_F)^T \Delta\mathbf{u}^{(0)} = 0
$$

Solving:

$$
\delta\lambda = -\frac{\delta\mathbf{u}_R^T \Delta\mathbf{u}^{(0)}}{\delta\mathbf{u}_F^T \Delta\mathbf{u}^{(0)}} = \frac{\delta\mathbf{q}^T \Delta\mathbf{u}}{\delta\mathbf{f}^T \Delta\mathbf{u}}
$$

This is computationally efficient and robust for most structural problems.

---

## 5. Tangent Stiffness Matrix

### 5.1 Material tangent

From the constitutive relation $\boldsymbol{\sigma} = \boldsymbol{\sigma}(\boldsymbol{\varepsilon})$:

$$
\mathbf{K}_T^{mat} = \int_V \mathbf{B}^T \mathbf{D}_T \mathbf{B} \, dV
$$

where $\mathbf{D}_T = \frac{\partial \boldsymbol{\sigma}}{\partial \boldsymbol{\varepsilon}}$ is the tangent modulus.

For elasticity: $\mathbf{D}_T = \mathbf{D}$ (constant)

For plasticity: $\mathbf{D}_T = \mathbf{D}^{ep}$ (consistent elastoplastic tangent)

### 5.2 Geometric stiffness

For large deformations with the Green-Lagrange strain:

$$
\mathbf{K}_T = \mathbf{K}_T^{mat} + \mathbf{K}_G
$$

The geometric stiffness arises from the variation of the strain-displacement matrix:

$$
\mathbf{K}_G = \int_V \mathbf{G}^T \boldsymbol{\tau} \mathbf{G} \, dV
$$

where $\mathbf{G}$ contains shape function derivatives and $\boldsymbol{\tau}$ is the stress matrix.

---

## 6. Convergence Criteria

### 6.1 Residual norm

Force equilibrium:

$$
\frac{\|\mathbf{R}^{(k)}\|}{\|\mathbf{F}^{ext}\|} < \epsilon_R
$$

### 6.2 Displacement increment norm

Solution stationarity:

$$
\frac{\|\Delta\mathbf{u}^{(k)}\|}{\|\mathbf{u}^{(k)}\|} < \epsilon_u
$$

### 6.3 Energy norm

Combined criterion:

$$
\frac{|\Delta\mathbf{u}^{(k)} \cdot \mathbf{R}^{(k)}|}{|\Delta\mathbf{u}^{(0)} \cdot \mathbf{R}^{(0)}|} < \epsilon_E
$$

The `femlabpy` helper `rnorm()` evaluates the residual on constrained DOFs, matching legacy conventions.

---

## 7. Implementation in femlabpy

### 7.1 Orthogonal residual solver

`solve_nlbar(X, T, G, C, P, ...)` implements the nonlinear bar solver:

| Step | Action |
|------|--------|
| 1 | Initialize $\mathbf{u} = \mathbf{0}$, $\mathbf{f} = \mathbf{0}$ |
| 2 | Build tangent with `kbar` |
| 3 | Apply BCs with `setbc` |
| 4 | Compute reference increment $\Delta\mathbf{u}_0$ |
| 5 | Orthogonal residual correction loop |
| 6 | Commit converged state |
| 7 | Recover reactions with `reaction` |

Returns: `u`, `q`, `S`, `E`, `R`, `f`, `U_path`, `F_path`

### 7.2 Elastoplastic solver

`solve_plastic(X, T, G, C, P, ...)` handles Q4 elastoplastic problems:

**Geometry options:**
- Plane stress: `kq4eps`, `qq4eps`
- Plane strain: `kq4epe`, `qq4epe`

**Material options:**
- `material_type=1`: von Mises
- `material_type=2`: Drucker-Prager

**State arrays:**
- `u`, `du`: displacement history
- `f`, `df`: external force state
- `S`, `E`: committed element stress/strain

### 7.3 Code reading order

```
solve_nlbar()
    └── kbar, qbar (element routines)
    └── setbc, setload (constraints)
    └── reaction (post-processing)

solve_plastic()
    └── kq4eps/kq4epe, qq4eps/qq4epe
    └── _solve_plastic_system()
```
