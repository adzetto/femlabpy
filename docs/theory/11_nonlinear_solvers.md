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

# Chapter 11: Nonlinear Solvers

This chapter covers the two nonlinear control flows implemented in the codebase:
the orthogonal-residual load-stepping path for nonlinear bars and the Newton
style elastoplastic solve for Q4 elements. Both routines are legacy-compatible
wrappers around the actual element kernels and constraint helpers.

## Orthogonal Residual Solver

`solve_nlbar(X, T, G, C, P, ...)` reproduces the nonlinear bar workflow used by
the old FemLab classroom examples. The function is not a generic arc-length
framework. It is a concrete solver with a specific data flow:

1. start from zero displacement, zero internal force, and the applied nodal load
   table;
2. build a tangent stiffness matrix with `kbar`;
3. apply displacement constraints with `setbc`;
4. solve a reference increment `du0`;
5. iterate with an orthogonal residual correction until convergence;
6. commit the load and displacement path;
7. recover reactions with `reaction`.

The returned dictionary contains the converged fields and the load path:

- `u`
- `q`
- `S`
- `E`
- `R`
- `f`
- `U_path`
- `F_path`

### What the loop actually does

The solver keeps three vectors in play:

- `u`, the accumulated displacement state,
- `du`, the current incremental step,
- `f`, the accumulated external force state.

For each load step it computes a fresh tangent matrix, solves a reference
increment, then corrects that increment using the residual from
`qbar(..., u + du)`. The ratio

`((dq.T @ du) / (df.T @ du))`

is the scalar update that moves the step along the orthogonal-residual path.
If the corrector does not converge before `i_max`, the solver restarts with a
smaller step.

### Why the residual norm matters

The helper `rnorm()` measures the residual on the constrained system in the same
legacy convention as the old code. That means the convergence test is tied to the
active degrees of freedom rather than a raw global vector norm.

### Practical meaning of the inputs

- `X` is the node coordinate table.
- `T` is the bar topology table.
- `G` stores material and section data.
- `C` stores prescribed displacements.
- `P` stores nodal loads.
- `plotdof` is one-based, because the original benchmark scripts were one-based.

## Elastoplastic Solver

`solve_plastic(X, T, G, C, P, ...)` handles the legacy Q4 elastoplastic cases.
The code keeps the constitutive update inside the element routines and uses the
solver only for load stepping and equilibrium enforcement.

The solver supports two geometric branches:

- plane stress through `kq4eps` and `qq4eps`;
- plane strain through `kq4epe` and `qq4epe`.

It also supports two material laws selected by `material_type`:

- `1` for von Mises;
- `2` for Drucker-Prager.

### Plane strain and plane stress

The `_solve_plastic_system()` helper chooses a symmetry-aware dense fallback for
plane strain and the generic linear solver otherwise. That is a compatibility
choice, not a new constitutive model.

### State flow

The plastic driver keeps the following state arrays:

- `u` and `du` for the displacement history,
- `f` and `df` for the external force state,
- `S` and `E` for the committed element stress and strain-like history values.

Inside each corrector iteration it:

1. rebuilds the tangent matrix from the current state;
2. applies the boundary conditions;
3. computes the internal force vector with the Q4 stress update routine;
4. forms the residual `dq = q - f`;
5. updates the step with the orthogonal residual correction;
6. commits the element state when the step converges.

The returned dictionary mirrors the bar solver and includes the same path
vectors together with the element state arrays.

### What to read in the source

If you want to modify this solver, read the code in this order:

1. `solve_nlbar()`
2. `solve_plastic()`
3. `_solve_plastic_system()`
4. `kbar`, `qbar`, `kq4epe`, `kq4eps`, `qq4epe`, `qq4eps`
5. `setbc()`, `setload()`, and `reaction()`

That order matches the actual control flow and makes the restart behavior much
easier to follow.

