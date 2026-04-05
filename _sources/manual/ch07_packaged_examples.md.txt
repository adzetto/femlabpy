# Chapter 7: Packaged Workflows and Examples

This chapter mirrors the packaged drivers under `src/femlabpy/examples/`. The goal is not to hide the mechanics behind a one-line wrapper, but to show how the data dictionaries, assemblers, solvers, and post-processing steps fit together.

The examples below are intentionally practical:

- the nonlinear truss example shows how the legacy orthogonal-residual driver walks a snap-through path,
- the wave example shows why explicit dynamics wants a lumped mass matrix and a stable time step,
- the imported-mesh example shows how Gmsh data, direct constraints, and reaction recovery connect in one workflow.

## 7.1 Geometrically Nonlinear Truss Snap-Through

The `bar01` benchmark is the smallest example in the repository that still behaves like a real nonlinear structural problem. The geometry is shallow enough that the response does not stay on a single equilibrium branch, so a plain load-controlled solve is not robust. That is why the example uses `solve_nlbar`, which follows the legacy orthogonal-residual continuation strategy.

The packaged data come from `femlabpy.examples.bar01_data()`. The fields matter:

- `X` stores node coordinates,
- `T` stores bar connectivity and material ids,
- `G` stores the bar material table,
- `C` stores prescribed displacements,
- `P` stores the reference nodal load vector,
- `no_loadsteps`, `i_max`, `i_d`, and `TOL` control the continuation loop.

```python
import femlabpy as fp
from femlabpy.examples import bar01_data

data = bar01_data()

result = fp.solve_nlbar(
    data["X"],
    data["T"],
    data["G"],
    data["C"],
    data["P"],
    no_loadsteps=int(data["no_loadsteps"][0, 0]),
    i_max=int(data["i_max"][0, 0]),
    i_d=int(data["i_d"][0, 0]),
    tol=float(data["TOL"][0, 0]),
    plotdof=int(data["plotdof"][0, 0]),
)

u_path = result["U_path"][:, 0]
f_path = result["F_path"][:, 0]
```

What this code is doing is straightforward once you align it with the solver:

- `bar01_data()` gives you the exact legacy deck that the original MATLAB example used.
- `solve_nlbar(...)` assembles the tangent stiffness with `kbar`, evaluates internal forces with `qbar`, and applies displacement constraints with `setbc`.
- Each load step starts from a predictor displacement increment, then corrects that increment until the residual norm is small enough.
- `U_path` and `F_path` record the monitored degree of freedom, which is what you plot to see the snap-through curve.

The important point is that the solver is not trying to force the structure through a bad equilibrium path. It keeps relinearizing the tangent problem and lets the load factor change sign when the response demands it. That is exactly what a snap-through benchmark needs.

If you prefer the packaged driver instead of the lower-level call, `femlabpy.examples.run_bar01_nlbar(plot=True)` returns the same result plus ready-made figures.

## 7.2 Explicit Bar Wave Propagation

The wave example exists to show a different kind of numerical behavior. Here the goal is not geometric nonlinearity, but a time discretization that can follow a sharp wave front without smearing it away. That is why the example uses the explicit central-difference solver together with a lumped mass matrix.

The packaged driver is `femlabpy.examples.run_dynamic_wave()`. Its output includes the exact wave speed from the material data, the critical time step, the time step actually used, and a `TimeHistory` object with the full displacement history.

```python
from femlabpy.examples import run_dynamic_wave

result = run_dynamic_wave(dt_factor=0.8, plot=True)

print("critical dt:", result["dt_cr"])
print("used dt:", result["dt_used"])
print("exact wave speed:", result["c_wave_exact"])
print("measured wave speed:", result["c_wave_measured"])
```

The internal workflow is the part worth understanding:

- `dynamic_wave_data()` builds a 1D bar mesh, the material table, and the reference wave speed `c_wave = sqrt(E / rho)`.
- `_assemble_bar_system()` loops over the elements and accumulates the bar stiffness and lumped mass into the global matrices.
- `critical_timestep(K, M)` estimates the stability limit for the explicit method, then the example multiplies that limit by `dt_factor`.
- `solve_central_diff(...)` advances the response using the diagonal mass entries directly, so there is no matrix factorization in the time loop.
- the pulse is injected through a small local force function `p_func(t)`, which keeps the load history explicit and easy to inspect.

The history arrays are stored in time-major form. `TimeHistory.u` is organized so `result.u[i, j]` means "displacement at step `i`, degree of freedom `j`". That makes it easy to plot a waterfall snapshot, a midpoint trace, or a response history at a single node.

The example also checks the measured arrival time against the analytical estimate `L / c_wave`. That is useful because it turns the script into a validation case instead of just a pretty animation.

## 7.3 Imported Mesh With Direct Constraints

The triangular Gmsh example shows the most typical non-legacy workflow in the repository: read a mesh, assemble a system, apply prescribed displacements, solve, then recover reactions and stresses.

```python
from femlabpy.examples import gmsh_triangle_data
from femlabpy import init, kt3e, qt3e, reaction, setload, solve_lag

data = gmsh_triangle_data()

K, p, q = init(data["X"].shape[0], data["dof"], use_sparse=False)
K = kt3e(K, data["T"], data["X"], data["G"])
p = setload(p, data["P"])
u = solve_lag(K, p, data["C"], data["dof"])
q, S, E = qt3e(q, data["T"], data["X"], data["G"], u)
R = reaction(q, data["C"], data["dof"])
```

This example is useful because it makes the data flow explicit:

- `gmsh_triangle_data()` uses `load_gmsh` under the hood and returns a normal FemLab-style dictionary.
- `init(...)` creates the global stiffness, load, and internal-force arrays with the right size for the mesh.
- `kt3e(...)` scatters the element stiffness contributions into the global matrix.
- `setload(...)` writes the nodal loads into the global right-hand side vector.
- `solve_lag(...)` enforces the prescribed displacements by augmenting the system with Lagrange multipliers.
- `qt3e(...)` recomputes the internal forces and the element stress/strain fields from the solved displacement vector.
- `reaction(...)` extracts the support reactions from the internal-force vector so you can check equilibrium.

This is the pattern to follow whenever you want to work from an imported mesh instead of one of the packaged benchmark decks. The solver remains transparent: every array is still available, and every step is still inspectable.

## 7.4 What To Read Next

The next step is to look at the source code for the exact functions these examples call:

- `src/femlabpy/solvers.py` for `solve_nlbar`,
- `src/femlabpy/dynamics.py` for `solve_central_diff`, `critical_timestep`, and the load builders,
- `src/femlabpy/examples/legacy_cases.py` for the packaged nonlinear benchmark decks,
- `src/femlabpy/examples/dynamic_wave.py` for the wave-propagation validation script,
- `src/femlabpy/examples/gmsh_triangle.py` for the imported-mesh workflow.

If you read the examples and the source together, the library stops looking like a black box. You can see exactly where the input tables go, how the global matrices are assembled, and how the solver returns the response history you plot afterwards.
