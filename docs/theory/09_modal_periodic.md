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

# Modal analysis and periodic boundaries

This chapter covers two workflows that sit on top of the standard assembly and
solver path:

- extracting modal data from the generalized eigenvalue problem,
- enforcing periodic constraints for representative volume element (RVE)
  calculations and homogenization.

The code paths are separate, but they share the same global matrix conventions
and one-based legacy tables. `solve_modal()` works with reduced free-DOF blocks;
the periodic module builds a constraint system and solves the resulting
saddle-point problem with Lagrange multipliers.

## 1. Modal Analysis

### 1.1 Generalized eigenvalue problem

The free vibration equation without damping is:

$$
\mathbf{M}\ddot{\mathbf{u}} + \mathbf{K}\mathbf{u} = \mathbf{0}
$$

Assuming harmonic motion $\mathbf{u}(t) = \boldsymbol{\phi} e^{i\omega t}$:

$$
-\omega^2 \mathbf{M}\boldsymbol{\phi} e^{i\omega t} + \mathbf{K}\boldsymbol{\phi} e^{i\omega t} = \mathbf{0}
$$

This yields the generalized eigenvalue problem:

$$
\mathbf{K}\boldsymbol{\phi} = \omega^2 \mathbf{M}\boldsymbol{\phi}
$$

where $\omega$ is the circular natural frequency and $\boldsymbol{\phi}$ is the mode shape.

### 1.2 Orthogonality of modes

For distinct eigenvalues $\omega_m^2 \neq \omega_n^2$, the mode shapes satisfy:

$$
\boldsymbol{\phi}_m^T \mathbf{M} \boldsymbol{\phi}_n = 0, \quad \boldsymbol{\phi}_m^T \mathbf{K} \boldsymbol{\phi}_n = 0 \quad (m \neq n)
$$

**Proof:** From the eigenvalue equations:
$$
\mathbf{K}\boldsymbol{\phi}_m = \omega_m^2 \mathbf{M}\boldsymbol{\phi}_m, \quad \mathbf{K}\boldsymbol{\phi}_n = \omega_n^2 \mathbf{M}\boldsymbol{\phi}_n
$$

Pre-multiply the first by $\boldsymbol{\phi}_n^T$ and the second by $\boldsymbol{\phi}_m^T$:
$$
\boldsymbol{\phi}_n^T \mathbf{K}\boldsymbol{\phi}_m = \omega_m^2 \boldsymbol{\phi}_n^T \mathbf{M}\boldsymbol{\phi}_m
$$
$$
\boldsymbol{\phi}_m^T \mathbf{K}\boldsymbol{\phi}_n = \omega_n^2 \boldsymbol{\phi}_m^T \mathbf{M}\boldsymbol{\phi}_n
$$

Since $\mathbf{K}$ and $\mathbf{M}$ are symmetric, transposing gives equal LHS. Subtracting:
$$
(\omega_m^2 - \omega_n^2) \boldsymbol{\phi}_n^T \mathbf{M}\boldsymbol{\phi}_m = 0
$$

For $\omega_m \neq \omega_n$: $\boldsymbol{\phi}_n^T \mathbf{M}\boldsymbol{\phi}_m = 0$ ∎

### 1.3 Mass normalization

Modes are typically normalized such that:

$$
\boldsymbol{\phi}_n^T \mathbf{M} \boldsymbol{\phi}_n = 1
$$

This implies:
$$
\boldsymbol{\phi}_n^T \mathbf{K} \boldsymbol{\phi}_n = \omega_n^2
$$

### 1.4 Modal participation factors

The participation factor for mode $n$ in direction $j$ quantifies how much the mode contributes to the response under uniform excitation:

$$
\Gamma_{n,j} = \frac{\boldsymbol{\phi}_n^T \mathbf{M} \mathbf{r}_j}{\boldsymbol{\phi}_n^T \mathbf{M}\boldsymbol{\phi}_n}
$$

where $\mathbf{r}_j$ is a unit influence vector (ones for DOFs in direction $j$, zeros elsewhere).

For mass-normalized modes ($\boldsymbol{\phi}_n^T \mathbf{M}\boldsymbol{\phi}_n = 1$):

$$
\Gamma_{n,j} = \boldsymbol{\phi}_n^T \mathbf{M} \mathbf{r}_j
$$

### 1.5 Effective modal mass

The effective modal mass in direction $j$ is:

$$
M_{eff,n,j} = \Gamma_{n,j}^2 \cdot (\boldsymbol{\phi}_n^T \mathbf{M}\boldsymbol{\phi}_n) = \Gamma_{n,j}^2
$$

The sum of effective masses equals the total mass:

$$
\sum_{n=1}^{N} M_{eff,n,j} = \mathbf{r}_j^T \mathbf{M} \mathbf{r}_j = M_{total,j}
$$

### 1.6 Implementation details

`solve_modal(K, M, n_modes=10, C_bc=None, dof=2, sigma=0.0)` returns a `ModalResult` with:

- `eigenvalues`: $\omega^2$
- `omega`: circular frequencies (rad/s)
- `freq_hz`: $f = \omega/(2\pi)$ (Hz)
- `period`: $T = 1/f$ (seconds)
- `mode_shapes`: mass-normalized eigenvectors
- `participation`: $\Gamma_{n,j}$ per direction
- `effective_mass`: $M_{eff,n,j}$ per direction

The solver reduces the system by removing constrained DOFs before solving, then expands results back to full size.
it.

## 2. Periodic boundaries

The periodic module is organized around three steps:

1. identify matching nodes on opposite faces,
2. build the linear constraint system,
3. solve the augmented problem and average the result.

### Pairing nodes

`find_periodic_pairs(X, axis, tol=1e-6)` compares the minimum and maximum faces
of a domain along a selected axis and returns 1-based node pairs. The
implementation uses `cKDTree` when available and the node count is large enough;
otherwise it falls back to a direct distance search. The tolerance is scaled by
the domain size, which makes the matching rule stable across different units.

`find_all_periodic_pairs(X, periodic_axes, tol=1e-6)` repeats that process for
multiple axes and returns a dictionary keyed by axis.

`check_periodic_mesh(X, axis, tol=1e-6)` returns a plain report dictionary with
`valid`, `n_left`, `n_right`, `max_mismatch`, and `message`. This is the first
function to call if a mesh should be periodic but pairing fails.

### Fixing the rigid body mode

Pure periodic constraints do not remove rigid translations by themselves.
`fix_corner(X, C_existing, dof)` adds zero-displacement rows at the corner node
closest to the minimum coordinate corner so a separate direct solve does not
remain singular. `solve_periodic()` does not take a boundary-condition table;
it handles the periodic augmentation through `solve_lag_general()`.

### Constraint matrices and macro strain

`periodic_constraints(X, pairs, dof, eps_macro=None)` converts the node pairs
into the algebraic system

$$
\mathbf{G}\mathbf{u} = \mathbf{Q}
$$

where each row enforces one DOF pair. The right node gets `+1`, the left node
gets `-1`, and the macro strain contribution is written into `Q`. The Voigt
strain vector is converted to a tensor by `_voigt_to_tensor()`, so the shear
components are split consistently between the off-diagonal tensor entries.

`apply_macro_strain()` is just a convenience wrapper that returns the `Q`
vector from the same logic.

`solve_periodic(K, p, X, pairs, dof, eps_macro=None, return_lagrange=False)`
calls `periodic_constraints()` and then passes `K`, `p`, `G`, and `Q` to
`solve_lag_general()`. That means the periodic solve is not a special-case
solver. It is a standard saddle-point solve with a different constraint matrix.

## 3. Volume averages and homogenization

`volume_average_stress()` and `volume_average_strain()` are the reduction layer
that converts element-level fields into domain averages. The current
implementation supports `t3` and `q4` elements.

- For `t3`, the code reuses the triangular batch geometry helper and weights
  element values by area.
- For `q4`, the code averages the four Gauss-point values and estimates the
  area with a shoelace formula on the nodal polygon.

`homogenize(K, T, X, G_mat, pairs, dof, element_type="q4")` then drives the
RVE through unit macro-strain cases and assembles the effective constitutive
matrix column by column. In 2D the result is `3 x 3`; in 3D it is `6 x 6`.

The important part of the implementation is the flow:

1. build a unit macro strain vector,
2. call `solve_periodic()` with zero external load,
3. average the resulting stress with `volume_average_stress()`,
4. insert the result into one column of `C_eff`.

That structure is why the function is easy to extend and easy to audit.

## 4. Practical workflow

The normal sequence for an RVE solve is:

```python
import numpy as np
from femlabpy.periodic import (
    check_periodic_mesh,
    find_periodic_pairs,
    solve_periodic,
    homogenize,
)

report = check_periodic_mesh(X, axis=0)
if not report["valid"]:
    raise ValueError(report["message"])

pairs_x = find_periodic_pairs(X, axis=0)
pairs_y = find_periodic_pairs(X, axis=1)
pairs = np.vstack([pairs_x, pairs_y])
u = solve_periodic(K, p_zero, X, pairs, dof=2, eps_macro=np.array([1.0, 0.0, 0.0]))
C_eff = homogenize(K, T, X, G_mat, pairs, dof=2, element_type="q4")
```

The same pattern is used in the repository examples: verify the mesh, pair the
faces, solve the constrained system, and average the stress response into an
effective stiffness tensor. If you need a separate reference constraint for a
non-periodic solve, add `fix_corner()` to your own `C` table before calling the
linear solver.
