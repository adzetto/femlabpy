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

## 1. Modal analysis

The modal solver works on the generalized eigenvalue problem

$$
\mathbf{K}\boldsymbol{\phi} = \omega^2 \mathbf{M}\boldsymbol{\phi}
$$

where `omega` is the circular natural frequency and `phi` is the mode shape.
`solve_modal(K, M, n_modes=10, C_bc=None, dof=2, sigma=0.0)` returns a
`ModalResult` dataclass with these fields:

- `eigenvalues`: the squared frequencies, `omega^2`,
- `omega`: circular frequencies in rad/s,
- `freq_hz`: frequencies in Hz,
- `period`: modal periods in seconds,
- `mode_shapes`: full-size mode shape vectors with zeros at constrained DOFs,
- `participation`: modal participation factors per direction,
- `effective_mass`: effective modal mass per direction.

### How the solver reduces the system

The implementation does not solve the full constrained system directly. It
first builds a boolean mask of free DOFs from `C_bc` with `_get_free_dofs()`,
then extracts the reduced matrices with `_reduce_system()`. The reduction uses
`np.ix_` so both the stiffness and mass matrices are sliced consistently:

```python
free_idx = np.where(free_mask)[0]
K_red = K[np.ix_(free_idx, free_idx)]
M_red = M[np.ix_(free_idx, free_idx)]
```

That matters because the solver is preserving the generalized eigenvalue
structure, not just deleting rows and columns arbitrarily.

For small systems, `solve_modal()` uses dense `scipy.linalg.eigh`. For larger
systems it switches to `scipy.sparse.linalg.eigsh` with shift-invert. After the
eigensolve, the eigenvectors are mass-normalized on the reduced system and then
expanded back to the full DOF space.

### Participation factors and effective mass

The internal helper `_modal_participation()` builds a unit influence vector for
each spatial direction and evaluates

$$
\Gamma_{n,j} = \frac{\phi_n^T \mathbf{M} r_j}{\phi_n^T \mathbf{M}\phi_n}
$$

For mass-normalized modes, the denominator is one, so the effective modal mass
reduces to the square of the participation factor. The code uses the full mass
matrix after expansion so the returned arrays are aligned with the original DOF
numbering.

### Plotting modes

`plot_modes(T, X, phi, dof, mode_indices=None, scale=1.0)` deforms the mesh for
one or more modes and is intentionally separate from the solver. The solver is
responsible for modal data; the plotting helper is responsible for presenting
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
