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

# Structural Dynamics and Time Integration

This chapter explains the solver side of `femlabpy`: mass and damping models,
load histories, time integration, frequency response functions, and the
post-processing helpers that operate on the resulting time histories.

The semidiscrete equation of motion is

$$
\mathbf{M}\ddot{\mathbf{u}}(t) + \mathbf{C}\dot{\mathbf{u}}(t) + \mathbf{K}\mathbf{u}(t) = \mathbf{p}(t)
$$

where `M`, `C`, and `K` are the global mass, damping, and stiffness matrices,
and `p(t)` is a load vector returned by one of the load callables in
`femlabpy.dynamics`. The time-history solvers all return the same container,
`TimeHistory`, with rows ordered by time:

- `t`: `(nsteps + 1,)`
- `u`: `(nsteps + 1, ndof)`
- `v`: `(nsteps + 1, ndof)`
- `a`: `(nsteps + 1, ndof)`
- `energy`: `None` or a dictionary with `kinetic`, `strain`, and `total`

## 1. Mass and damping models

The mass matrix controls inertia, and the damping matrix controls dissipation.
The implementation keeps both concepts close to the global algebraic system so
the same routines can be used for dense and sparse workflows.

### Consistent and lumped mass

`femlabpy` uses the two mass patterns that are standard in structural
dynamics:

- consistent mass, which follows the same interpolation as the stiffness
  kernels;
- lumped mass, which compresses the row sums onto the diagonal.

The distinction matters because `solve_central_diff()` requires a diagonal mass
matrix. The helper `_get_lumped_diagonal()` in `src/femlabpy/dynamics.py`
checks this explicitly and raises an error if the matrix is not diagonal.

For stability estimates, `critical_timestep()` uses a power iteration on
`M^-1 K` and returns approximately `2 / omega_max`. That is the practical
quantity you want when deciding whether an explicit scheme is safe.

### Rayleigh and modal damping

`rayleigh_coefficients(omega1, omega2, zeta1, zeta2)` solves the 2-by-2 system
that matches two target modal damping ratios. The coefficients are then fed to
`rayleigh_damping(M, K, alpha, beta)`, which returns

$$
\mathbf{C} = \alpha \mathbf{M} + \beta \mathbf{K}
$$

If you want a damping matrix that exactly matches a selected set of modal
damping ratios, `modal_damping(M, omega, phi, zeta)` builds the dense Caughey
form

$$
\mathbf{C} = \mathbf{M}\Phi \operatorname{diag}(2\zeta_i\omega_i)\Phi^T\mathbf{M}
$$

using mass-normalized mode shapes.

## 2. Load histories

The load builders all return callables with the same shape convention:
`p(t) -> (ndof, 1)`.

- `constant_load(P)` returns a fixed vector.
- `ramp_load(P, t_ramp)` linearly ramps until the target time.
- `harmonic_load(P, omega, phase=0.0)` evaluates a sinusoid.
- `pulse_load(P, t_start, t_duration)` returns a rectangular pulse.
- `tabulated_load(P, time_table, value_table)` interpolates a scalar history
  with `np.interp`.
- `seismic_load(M, direction, accel_record, dt_record)` precomputes `-M @ r`
  once, then interpolates the ground acceleration record at runtime.

The seismic loader is worth calling out because it avoids repeated matrix-vector
products inside the time loop. The returned callable only interpolates the
record and scales the precomputed influence vector.

## 3. Newmark beta

`NewmarkParams` is a small convenience container for the integration constants.
The named constructors document the standard choices:

- `average_acceleration()` for the unconditionally stable trapezoidal rule,
- `linear_acceleration()` for the classical second-order variant,
- `central_difference()` for the explicit limiting case,
- `fox_goodwin()` for the higher-order SDOF variant.

`solve_newmark()` is the main implicit integrator. The implementation follows a
fixed sequence:

1. reshape `u0` and `v0` to column vectors,
2. determine constrained DOFs from `C_bc`,
3. build the effective stiffness matrix `K_eff = K + a0*M + a1*C`,
4. factorize `K_eff` once when SciPy is available,
5. compute the initial acceleration from `p(0) - C v0 - K u0`,
6. march forward in time, updating `u`, `v`, and `a`,
7. optionally compute kinetic and strain energy histories.

The solver writes history arrays in time-major order, so `result.u[n, i]` is the
displacement of DOF `i` at time step `n`. That is the layout used by the
plotting helpers and by downstream post-processing.

When `compute_energy=True`, the returned `TimeHistory.energy` dictionary stores
only the keys that the implementation currently computes: `kinetic`, `strain`,
and `total`. `plot_energy()` expects exactly that layout.

`solve_newmark_nl()` follows the same time-history pattern, but the update is
wrapped in a Newton iteration. It accepts a tangent stiffness callback and an
internal-force callback, then iterates until the residual satisfies the
specified tolerance or `max_iter` is reached.

## 4. Central difference and HHT alpha

`solve_central_diff()` is the explicit companion to Newmark beta. It requires a
lumped mass matrix, reconstructs the first step through backward extrapolation,
and then advances using the central difference formula. If a damping matrix is
present, the solver only uses its diagonal contribution.

This is the routine that pairs naturally with `critical_timestep()`: if the mass
is not diagonal or the time step is too large, the explicit scheme is not a good
fit.

`solve_hht()` implements the Hilber-Hughes-Taylor method. The parameter `alpha`
must stay in `[-1/3, 0]`; the solver derives

$$
\beta = \frac{(1 - \alpha)^2}{4}, \qquad \gamma = \frac{1}{2} - \alpha
$$

and then forms the modified effective stiffness and load. Compared with Newmark
average acceleration, the HHT method damps high-frequency content more
aggressively while keeping second-order accuracy for the low-frequency response.

## 5. Frequency response and plots

`compute_frf(M, C, K, input_dof, output_dof, freq_range, n_points=500)` sweeps a
frequency band in Hertz, converts each sample to angular frequency, and solves

$$
\left(\mathbf{K} - \omega^2\mathbf{M} + i\omega\mathbf{C}\right)\mathbf{u} = \mathbf{f}
$$

for a unit force applied at `input_dof`. The function returns the frequency grid
and the complex scalar transfer function at `output_dof`.

`plot_frf()` renders magnitude and phase in a two-panel figure and can mark
resonance peaks. `plot_time_history()` plots displacement, velocity, or
acceleration history for one or more DOFs. `plot_energy()` draws the current
`energy` dictionary and raises a clear error if the solver was not run with
energy tracking enabled.

## 6. Practical use

The standard workflow is short:

```python
from femlabpy.dynamics import (
    rayleigh_coefficients,
    rayleigh_damping,
    solve_newmark,
    plot_time_history,
)

alpha, beta = rayleigh_coefficients(omega1, omega2, 0.05, 0.05)
C = rayleigh_damping(M, K, alpha, beta)
history = solve_newmark(M, C, K, p_func, u0, v0, dt=0.01, nsteps=1000)
plot_time_history(history, dof_index=0, quantity="displacement")
```

The important part is the data flow: the load builder returns `p(t)`, the solver
turns it into a `TimeHistory`, and the plotting helpers read that object without
needing to know how the matrices were assembled.
