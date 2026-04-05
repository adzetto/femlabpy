# Chapter 9: Dynamic Response Workflows

This chapter focuses on the parts of `femlabpy` that are easiest to misuse if you treat them as black boxes: modal output, frequency response, seismic loading, time-history arrays, and the plotting helpers that sit on top of them. The examples below follow the actual implementation in `src/femlabpy/dynamics.py`, `src/femlabpy/modal.py`, `src/femlabpy/damping.py`, `src/femlabpy/plotting.py`, and `src/femlabpy/postprocess.py`.

## 9.1 Modal analysis workflow

Modal analysis in `femlabpy` solves the generalized eigenvalue problem

$$
\mathbf{K} \phi = \omega^2 \mathbf{M} \phi
$$

after constrained degrees of freedom have been removed from the system. The public result is a `ModalResult` object with these fields:

- `eigenvalues`: the raw eigenvalues, which are the squared circular frequencies.
- `omega`: circular frequencies in rad/s.
- `freq_hz`: frequencies in Hz.
- `period`: periods in seconds.
- `mode_shapes`: full-size mode shape vectors, padded back to the global DOF space with zeros at constrained DOFs.
- `participation`: modal participation factors per direction.
- `effective_mass`: effective modal masses per direction.

That output matters because the mode shapes returned by `solve_modal()` are already mapped back to the full system. You do not need to rebuild the constrained DOFs yourself unless you want a reduced representation for a custom postprocessor.

### How `solve_modal` reduces the system

The implementation in `src/femlabpy/modal.py` first converts the boundary-condition table into a boolean mask of free DOFs, then extracts the free-free blocks of `K` and `M` with `np.ix_`. This is why `np.ix_` matters here: it is doing row and column selection at the same time, so the solver only sees the dynamic part of the model.

```python
free_dofs = np.setdiff1d(np.arange(ndof), constrained_dofs)
K_free = K[np.ix_(free_dofs, free_dofs)]
M_free = M[np.ix_(free_dofs, free_dofs)]
```

For small systems the code uses dense `scipy.linalg.eigh`; for larger systems it switches to `scipy.sparse.linalg.eigsh`. After the solve, the eigenvectors are mass-normalized and expanded back into the full DOF space.

### What the modal output means in practice

The most useful fields are `freq_hz`, `period`, and `effective_mass`. `freq_hz` is what you normally report in a table. `period` is the same result expressed as `1 / f`. `effective_mass` tells you how much mass each mode mobilizes in each direction, which is the value you use when deciding whether enough modes were extracted for response-spectrum or time-history work.

`femlabpy.modal.plot_modes()` is the companion plotting helper. It takes the full mode-shape matrix, reshapes each mode by node and DOF, and overlays the deformed mesh on the original one.

### Example: cantilever modal study

```{code-block} python
import numpy as np
import matplotlib.pyplot as plt

import femlabpy as fp
from femlabpy.modal import solve_modal, plot_modes

# Build a small 2D cantilever model.
data = fp.canti()
T = data["T"]
X = data["X"]
G = data["G"]
C = data["C"]
dof = int(data["dof"])
nn = X.shape[0]
ndof = nn * dof

K = np.zeros((ndof, ndof), dtype=float)
M = np.zeros((ndof, ndof), dtype=float)
K = fp.kq4e(K, T, X, G)
M = fp.mq4e(M, T, X, G, lumped=False)

result = solve_modal(K, M, n_modes=5, C_bc=C, dof=dof)

print("mode  freq_hz  period  effective_mass_x")
for i in range(len(result.freq_hz)):
    print(
        f"{i + 1:>4}  {result.freq_hz[i]:>8.3f}  {result.period[i]:>7.4f}  "
        f"{result.effective_mass[i, 0]:>16.3f}"
    )

fig = plot_modes(T, X, result.mode_shapes, dof, mode_indices=[0, 1, 2], scale=0.25)
plt.show()
```

The two lines that matter most are the mass assembly and the modal solve. Everything else in this example is just there to show how the result object is used after the solve.

## 9.2 Frequency response workflow

`compute_frf()` evaluates the complex dynamic stiffness

$$
\mathbf{Z}(\omega) = \mathbf{K} - \omega^2 \mathbf{M} + i \omega \mathbf{C}
$$

at a set of frequencies and solves one linear system per frequency. The function returns two arrays:

- `freq_hz`: the sampled frequency axis in Hz.
- `H`: the complex transfer function between the selected input and output DOFs.

This is a scalar FRF, not a full matrix. The input and output DOF indices are global and 0-based.

`plot_frf()` then converts `H` into magnitude and phase plots. If `mark_peaks=True`, the helper uses `scipy.signal.find_peaks` to mark likely resonance peaks on the magnitude plot.

### Example: one input, one output

```{code-block} python
import numpy as np
import matplotlib.pyplot as plt

from femlabpy.dynamics import compute_frf, plot_frf

# Assume K, M, and C are already assembled and boundary conditions are applied.
freq_hz, H = compute_frf(
    M,
    C,
    K,
    input_dof=0,
    output_dof=0,
    freq_range=(0.1, 50.0),
    n_points=800,
)

fig = plot_frf(freq_hz, H, log_scale=True, mark_peaks=True)
plt.show()
```

If you only need the resonance locations, the sampled magnitude is usually enough. If you need a transfer function for later analysis, keep `H` and postprocess it yourself.

## 9.3 Seismic loading workflow

`seismic_load()` is a load-function factory. It does not solve the dynamic problem by itself. Instead, it precomputes the effective base-excitation vector

$$
\mathbf{p}(t) = -\mathbf{M} \mathbf{r} a_g(t)
$$

and returns a callable `p(t)` that interpolates the recorded ground acceleration at any requested time.

What matters in the code is that the function accepts either a dense or sparse mass matrix, multiplies it by the influence vector once, and then uses `np.interp` for the ground-motion history. That keeps the time integrator simple: the solver only needs to call `p_func(t_next)` at each step.

### Example: ground-motion input

```{code-block} python
import numpy as np
import matplotlib.pyplot as plt

from femlabpy.dynamics import seismic_load, solve_newmark, plot_time_history, plot_energy

dt_record = 0.01
n_points = 1000
t_record = np.arange(n_points) * dt_record
ag_g = 0.5 * np.sin(2.0 * np.pi * 2.5 * t_record) * np.exp(-0.2 * t_record)
ag = ag_g * 9.80665

# Reuse the same assembled model and boundary-condition table from the modal example.
C_damp = np.zeros_like(K)

# Influence vector: 1 for x-direction DOFs, 0 for all others.
influence = np.zeros(ndof, dtype=float)
influence[0::2] = 1.0

p_func = seismic_load(M, influence, ag, dt_record)

u0 = np.zeros((ndof, 1))
v0 = np.zeros((ndof, 1))

history = solve_newmark(
    M,
    C_damp,
    K,
    p_func,
    u0,
    v0,
    dt=dt_record,
    nsteps=n_points,
    C_bc=data["C"],
    dof=dof,
    compute_energy=True,
)
```

`solve_newmark()` returns a `TimeHistory` object. Its arrays are stored time-first:

- `history.t` has shape `(nsteps + 1,)`.
- `history.u`, `history.v`, and `history.a` have shape `(nsteps + 1, ndof)`.
- Row `0` is the initial state.
- Column `j` is the history of global DOF `j`.

That means a single DOF history is extracted as `history.u[:, target_dof]`, not by reshaping the array.

### Example: extracting one roof DOF

```{code-block} python
roof_node = np.argmax(X[:, 1])
roof_dof_x = roof_node * dof
u_roof = history.u[:, roof_dof_x]

plt.figure(figsize=(10, 4))
plt.plot(history.t, u_roof * 1000.0, color="tab:blue", linewidth=1.5)
plt.xlabel("Time (s)")
plt.ylabel("Displacement (mm)")
plt.title("Roof displacement history")
plt.grid(True, alpha=0.3)
plt.show()
```

`plot_time_history()` is the library helper for the same task. It accepts one DOF index or a list of indices and plots displacement, velocity, or acceleration directly from the `TimeHistory` object.

## 9.4 Energy and postprocessing workflow

When `compute_energy=True`, the solvers populate `history.energy` with a dictionary. In the current implementation the keys are:

- `kinetic`
- `strain`
- `total`

`plot_energy()` expects exactly that dictionary. If the solver was not run with energy tracking enabled, the helper raises a `ValueError` instead of guessing.

### Example: energy plot

```{code-block} python
fig, ax = plt.subplots(figsize=(10, 4))
plot_energy(history, ax=ax)
plt.show()
```

The implementation is intentionally simple. It does not attempt to reconstruct external work or damping work histories. If you need those, you should compute them explicitly from your loading history and state vectors.

### Example: using the built-in time-history plotter

```{code-block} python
fig, ax = plt.subplots(figsize=(10, 4))
plot_time_history(history, [roof_dof_x], quantity="displacement", ax=ax)
plt.show()
```

This helper is useful when you want the standard library behavior and do not want to write array slicing and axis labeling by hand.

## 9.5 Recommended workflow

A practical dynamic workflow in `femlabpy` usually looks like this:

1. Assemble `K`, `M`, and optionally `C`.
2. Run `solve_modal()` to check periods, mode shapes, and modal mass participation.
3. Use `compute_frf()` when you need steady-state response versus frequency.
4. Build a `p_func` with `constant_load()`, `harmonic_load()`, `tabulated_load()`, or `seismic_load()`.
5. Run `solve_newmark()`, `solve_hht()`, `solve_central_diff()`, or `solve_newmark_nl()` depending on the problem.
6. Postprocess with `plot_time_history()`, `plot_energy()`, and `plot_frf()`.

That sequence keeps the solver, the history arrays, and the postprocessing helpers aligned with the actual data structures used by the package.
