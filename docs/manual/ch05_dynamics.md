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

This chapter explains how `femlabpy` moves from static equilibrium to time-dependent response. The library follows the standard semi-discrete form

$$
\mathbf{M}\ddot{\mathbf{u}}(t) + \mathbf{C}\dot{\mathbf{u}}(t) + \mathbf{K}\mathbf{u}(t) = \mathbf{p}(t)
$$

but the important part for users is how the code is organized: loads are passed in as callables, constrained degrees of freedom are handled explicitly, and the solver results are stored in a `TimeHistory` container with arrays shaped as `(nsteps + 1, ndof)`.

The code in this chapter matches the actual `femlabpy.dynamics`, `femlabpy.damping`, and `femlabpy.modal` modules, so the narrative below is written around the real API rather than a textbook-only presentation.

## 5.1 Data layout and solver conventions

`femlabpy` uses the same indexing conventions across assembly, modal analysis, and time integration:

- Node numbering in user-facing tables is one-based.
- DOF indices returned by the dynamic solvers are zero-based.
- Constraint tables such as `C_bc` are interpreted with one-based node and local-DOF columns.
- Load functions return column vectors with shape `(ndof, 1)`.

The `TimeHistory` dataclass stores:

- `t`: time stamps with shape `(nsteps + 1,)`
- `u`: displacement history with shape `(nsteps + 1, ndof)`
- `v`: velocity history with shape `(nsteps + 1, ndof)`
- `a`: acceleration history with shape `(nsteps + 1, ndof)`
- `energy`: optional dictionary with the keys `kinetic`, `strain`, and `total`

That shape choice is deliberate. Each row is one time step, which makes slicing and plotting straightforward:

```python
roof_history = result.u[:, roof_dof]
```

## 5.2 Mass matrices

The mass matrix controls inertia. In practice, the choice is not only mathematical but also algorithmic.

### Consistent mass

The consistent mass matrix uses the same shape functions as the stiffness matrix:

$$
\mathbf{M}_e = \int_{\Omega_e} \rho \mathbf{N}^T \mathbf{N}\, d\Omega
$$

This formulation is the better choice for modal analysis and implicit transient analysis because it preserves coupling between DOFs. It is also what you want when the accuracy of the first few modes matters.

### Lumped mass

The lumped mass matrix is diagonal. In `femlabpy`, that matters for two reasons:

- `solve_central_diff` requires a diagonal mass matrix and rejects a non-diagonal input.
- The diagonal form makes explicit updates cheap because each DOF can be updated independently.

The code accepts either a diagonal vector or a full diagonal matrix. Internally it extracts the diagonal and checks that the off-diagonal norm is negligible before continuing.

For example, a lumped mass is the right choice when the model is used with very small time steps and the goal is speed rather than spectral fidelity.

## 5.3 Modal analysis and damping

### Modal analysis

`solve_modal(K, M, ...)` solves the generalized eigenvalue problem

$$
\mathbf{K}\boldsymbol{\phi} = \omega^2 \mathbf{M}\boldsymbol{\phi}
$$

but the implementation has a few practical steps that are easy to miss:

- Constrained DOFs are removed before the eigensolve.
- Small systems use dense `scipy.linalg.eigh`.
- Larger systems switch to `scipy.sparse.linalg.eigsh`.
- Mode shapes are mass-normalized after the eigensolve.
- The returned mode shapes are expanded back to the full DOF size with zeros at fixed DOFs.

That means the output is ready to use directly in plotting or in downstream response calculations.

The returned `ModalResult` also includes participation factors and effective modal mass. Those are computed from the full mass matrix and the full-size mode shapes, not from the reduced eigensystem. This is important because the influence vectors are built per active direction, and the modal mass calculation must respect the actual global numbering.

### Rayleigh damping

`rayleigh_coefficients(omega1, omega2, zeta1, zeta2)` solves the two-point calibration problem for

$$
\mathbf{C} = \alpha \mathbf{M} + \beta \mathbf{K}
$$

The solver is just a `2 x 2` linear system, but the implementation has a useful interpretation:

- `alpha` controls low-frequency damping.
- `beta` controls high-frequency damping.
- Frequencies between the calibration points get intermediate damping.

`rayleigh_damping(M, K, alpha, beta)` then assembles the matrix. It preserves sparsity when either input is sparse, which keeps the transient solvers practical on larger models.

### Modal damping

`modal_damping(M, omega, phi, zeta)` is a denser, more explicit alternative. It is useful for small reduced systems or validation cases where you want to prescribe damping mode by mode. The function constructs the damping matrix directly from mass-normalized modes, so it is not the first choice for large production models.

## 5.4 Load functions and frequency response

The transient solvers in `femlabpy` do not take a fixed load vector. They take a callable, which keeps the API simple and makes time interpolation local to the load generator.

### Load builders

`constant_load`, `ramp_load`, `harmonic_load`, `pulse_load`, and `tabulated_load` all return a function with the signature `p(t) -> (ndof, 1)` column vector.

This design avoids rebuilding the load array outside the time loop and makes the solver code read naturally:

```python
p_next = p_func(t_next)
```

`tabulated_load` uses `numpy.interp` so the input can be a sampled history. `seismic_load` does the same thing for base excitation, but first multiplies the mass matrix by the influence vector:

$$
\mathbf{p}_{eff}(t) = -\mathbf{M}\mathbf{r} a_g(t)
$$

This matters because `seismic_load` is not a generic force wrapper. It is specifically a base-excitation generator for relative-motion formulations.

### Frequency response function

`compute_frf(M, C, K, input_dof, output_dof, freq_range, n_points=500)` evaluates the complex dynamic stiffness

$$
\mathbf{Z}(\omega) = \mathbf{K} - \omega^2\mathbf{M} + i\omega\mathbf{C}
$$

and solves a dense linear system at each sampled frequency.

The practical consequences are:

- it is useful for small and medium systems;
- it returns the scalar transfer function between one input DOF and one output DOF;
- `plot_frf` then visualizes magnitude and phase, and can annotate peaks with `scipy.signal.find_peaks` when SciPy is available.

## 5.5 Newmark beta method

`solve_newmark` is the main implicit solver in the module. Its implementation is close to the textbook derivation, but a few code-level choices make it usable:

- The input load is a callable `p_func`.
- Initial acceleration is computed from equilibrium at `t = 0`.
- The effective stiffness matrix is formed once as `K + a0*M + a1*C`.
- If SciPy is available, the solver factorizes that matrix once with LU and reuses the factorization at every step.
- Boundary conditions are applied by zeroing the constrained rows and columns in the effective stiffness and by zeroing the constrained entries in the time histories.

The update formulas are the standard Newmark recurrences:

$$
\mathbf{u}_{n+1} = \mathbf{u}_n + \Delta t\,\mathbf{v}_n + \Delta t^2\left[\left(\frac{1}{2}-\beta\right)\mathbf{a}_n + \beta \mathbf{a}_{n+1}\right]
$$

$$
\mathbf{v}_{n+1} = \mathbf{v}_n + \Delta t\left[(1-\gamma)\mathbf{a}_n + \gamma \mathbf{a}_{n+1}\right]
$$

What the code adds on top of that is the efficient matrix handling:

```python
K_eff = K + a0 * M + a1 * C
K_factor = lu_factor(K_eff)
u_next = lu_solve(K_factor, p_eff)
```

That is the right implementation for a linear system with fixed `dt`, fixed matrices, and a load that changes only through `p_func(t)`.

### Parameter presets

`NewmarkParams` provides a few named presets that are worth documenting in the chapter:

- `average_acceleration()` corresponds to `beta = 0.25`, `gamma = 0.5`
- `linear_acceleration()` uses `beta = 1/6`, `gamma = 0.5`
- `central_difference()` encodes the explicit limit `beta = 0`, `gamma = 0.5`
- `fox_goodwin()` is included as a higher-order preset for undamped SDOF reference cases

These presets are useful as named intent, even when you still pass the scalar values into the solver.

### Energy tracking

When `compute_energy=True`, the solver records kinetic and strain energy at each step and stores them in `result.energy`. The current implementation does not try to reconstruct full work balance terms inside the solver; it keeps the energy output intentionally small and stable.

## 5.6 HHT alpha method

`solve_hht` is the same basic time-integration pipeline as Newmark, but with an extra high-frequency filter.

The code enforces `alpha` in the interval `[-1/3, 0]`, derives `beta` and `gamma` from it, and then modifies the effective stiffness and load so that the equilibrium equation is evaluated with the HHT weighting. In practice:

- `alpha = 0` reduces to standard Newmark average acceleration
- more negative values add more numerical dissipation at high frequencies
- the method is still second-order accurate for the low-frequency response

This is a good default when the model is stiff enough that spurious high-frequency oscillations are visible but you do not want to fully switch to an explicit formulation.

## 5.7 Central difference method

`solve_central_diff` is the explicit solver. Its implementation is stricter than the implicit solvers for a reason: it needs a diagonal mass matrix.

The code path does the following:

- extracts the lumped mass diagonal from a vector or diagonal matrix;
- checks that the matrix is actually diagonal;
- computes initial acceleration from the initial equilibrium state;
- forms the backward extrapolated displacement `u_{-1}`;
- steps forward using the explicit recurrence.

The main point is that the mass appears only as a diagonal scaling. That is why `solve_central_diff` can remain fast, but also why it is not a drop-in replacement for a consistent-mass transient analysis.

The stability limit is estimated with `critical_timestep(K, M)`, which uses power iteration on `M^{-1}K` to approximate the highest circular frequency and returns `2 / omega_max`.

## 5.8 Nonlinear Newmark

`solve_newmark_nl` extends the implicit Newmark pipeline with a Newton loop at every time step. Instead of assuming a fixed tangent stiffness, it accepts:

- `tangent_func(u, state)`
- `internal_force_func(u, state)`

That split is the right design for nonlinear materials or geometric nonlinearity because the element state can evolve during iteration.

The algorithm is:

1. Predict displacement and velocity.
2. Evaluate internal force and tangent stiffness at the trial state.
3. Form the residual `p - M a - C v - q`.
4. Solve the linear correction problem.
5. Update the trial state until the residual is below tolerance.

This is the part of the module to read if you want to understand how the library handles nonlinear dynamics without hiding the Newton iteration inside a black box.

## 5.9 Plotting and postprocessing

The plotting helpers are small but important because they match the solver output layout:

- `plot_time_history` plots displacement, velocity, or acceleration for selected DOFs
- `plot_energy` expects `compute_energy=True` and visualizes the energy channels in `TimeHistory.energy`
- `plot_modes` from `modal.py` plots deformed meshes from the modal result

These helpers are not just convenience wrappers. They encode the same data conventions as the solvers, so the examples in the manual can use them directly without reshaping the outputs by hand.

## 5.10 Example workflow

The following pattern is the one this repository is optimized for:

```python
from femlabpy.damping import rayleigh_coefficients, rayleigh_damping
from femlabpy.dynamics import seismic_load, solve_newmark
from femlabpy.modal import solve_modal

# Assemble K, M, and boundary conditions first.
result = solve_modal(K, M, n_modes=5, C_bc=C, dof=2)
alpha, beta = rayleigh_coefficients(result.omega[0], result.omega[1], 0.05, 0.05)
C_damp = rayleigh_damping(M, K, alpha, beta)

# Build a time-dependent load and solve the transient response.
p_func = seismic_load(M, influence_vector, accel_record, dt_record)
history = solve_newmark(
    M,
    C_damp,
    K,
    p_func,
    u0,
    v0,
    dt=dt_record,
    nsteps=len(accel_record) - 1,
    C_bc=C,
    dof=2,
)
```

That sequence mirrors the intended workflow in `femlabpy`: modal inspection first, damping calibration second, and transient analysis last.
