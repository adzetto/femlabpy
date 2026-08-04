# femlabpy

<!-- badges -->
![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white) ![TeX](https://img.shields.io/badge/TeX-008080?style=for-the-badge&logo=latex&logoColor=white) ![Scilab](https://img.shields.io/badge/Scilab-555555?style=for-the-badge) ![HTML](https://img.shields.io/badge/HTML-E34F26?style=for-the-badge&logo=html5&logoColor=white) ![MATLAB](https://img.shields.io/badge/MATLAB-0076A8?style=for-the-badge&logo=mathworks&logoColor=white) ![Julia](https://img.shields.io/badge/Julia-9558B2?style=for-the-badge&logo=julia&logoColor=white)

![NumPy](https://img.shields.io/badge/NumPy-013243?style=for-the-badge&logo=numpy&logoColor=white) ![SciPy](https://img.shields.io/badge/SciPy-8CAAE6?style=for-the-badge&logo=scipy&logoColor=white) ![Matplotlib](https://img.shields.io/badge/Matplotlib-11557C?style=for-the-badge) ![Express](https://img.shields.io/badge/Express-000000?style=for-the-badge&logo=express&logoColor=white) ![pytest](https://img.shields.io/badge/pytest-0A9EDC?style=for-the-badge&logo=pytest&logoColor=white) ![Qt](https://img.shields.io/badge/Qt-41CD52?style=for-the-badge&logo=qt&logoColor=white)

![last commit](https://img.shields.io/github/last-commit/adzetto/femlabpy?style=flat-square&color=informational) ![repo size](https://img.shields.io/github/repo-size/adzetto/femlabpy?style=flat-square&color=informational) ![top language](https://img.shields.io/github/languages/top/adzetto/femlabpy?style=flat-square) ![language count](https://img.shields.io/github/languages/count/adzetto/femlabpy?style=flat-square) ![license](https://img.shields.io/github/license/adzetto/femlabpy?style=flat-square&color=informational)


[![PyPI version](https://badge.fury.io/py/femlabpy.svg)](https://pypi.org/project/femlabpy/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/github/actions/workflow/status/adzetto/femlabpy/tests.yml?branch=main&label=tests)](https://github.com/adzetto/femlabpy/actions/workflows/tests.yml)
[![License](https://img.shields.io/github/license/adzetto/femlabpy)](LICENSE)

Python FEM library for teaching. Based on MATLAB/Scilab FemLab. Visit: https://adzetto.github.io/femlabpy

<br>
<p align="center">
  <b>
    <a href="#install">Install</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
    <a href="#complete-workflow-example">Workflow</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
    <a href="#quick-examples">Examples</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
    <a href="#api-reference">API Reference</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
    <a href="#dynamic-analysis--benchmarks">Benchmarks</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
    <a href="#development">Development</a>
  </b>
</p>
<br>

## Install

```bash
pip install femlabpy
```

Optional:
```bash
pip install "femlabpy[mesh]"  # Gmsh 4.x support
pip install "femlabpy[gui]"   # GUI tools
```

## Workflow Example

### Step 1: Create mesh with Gmsh

```python
import gmsh

gmsh.initialize()
gmsh.model.add("plate_with_hole")

# Create geometry
gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
gmsh.model.occ.addDisk(0.5, 0.5, 0, 0.1, 0.1)
gmsh.model.occ.cut([(2, 1)], [(2, 2)])
gmsh.model.occ.synchronize()

# Generate mesh
gmsh.model.mesh.generate(2)
gmsh.write("plate_hole.msh")
gmsh.finalize()
```

### Step 2: Load mesh into femlabpy

```python
import femlabpy as fp

mesh = fp.load_gmsh2("plate_hole.msh")
T = mesh.triangles  # element connectivity
X = mesh.positions  # node coordinates
```

### Step 3: Define material properties

```python
import numpy as np

E = 210000.0   # Young's modulus (MPa)
nu = 0.3       # Poisson's ratio
t = 1.0        # thickness (mm)

# G matrix: [E, nu, t, plane_stress_flag]
G = np.array([[E, nu, t, 1]])
```

### Step 4: Define boundary conditions

```python
# C matrix: [node, dof, value]
# Fix left edge (x=0): ux=0, uy=0
left_nodes = np.where(X[:, 0] < 0.01)[0] + 1
C = []
for n in left_nodes:
    C.append([n, 1, 0.0])  # ux = 0
    C.append([n, 2, 0.0])  # uy = 0
C = np.array(C)
```

### Step 5: Define loads

```python
# P matrix: [node, dof, value]
# Apply tension on right edge (x=1)
right_nodes = np.where(X[:, 0] > 0.99)[0] + 1
P = []
for n in right_nodes:
    P.append([n, 1, 100.0])  # Fx = 100 N
P = np.array(P)
```

### Step 6: Solve

```python
result = fp.elastic(T, X, G, C, P, dof=2)
u = result["u"]  # displacements
S = result["S"]  # stresses
```

### Step 7: Visualize

```python
fp.plotu(T, X, u, dof=2)
fp.plotelem(T, X)
```

## Quick Examples

### Cantilever beam (packaged example)

```python
import femlabpy as fp

data = fp.canti()
result = fp.elastic(data["T"], data["X"], data["G"], data["C"], data["P"], dof=2)
print("Max displacement:", result["u"].max())
print("Max stress:", result["S"].max())
```

### Potential flow

```python
# Q4 mesh
result_q4 = fp.flowq4(plot=True)
print("Temperature range:", result_q4["u"].min(), "to", result_q4["u"].max())

# T3 mesh
result_t3 = fp.flowt3(plot=True)
```

### Nonlinear truss (large deformation)

```python
from femlabpy.examples import run_bar01_nlbar

result = run_bar01_nlbar(plot=True)
print("Load path:", result["F_path"].ravel())
print("Displacement path:", result["U_path"].ravel())
```

### Plane stress plasticity (von Mises)

```python
from femlabpy.examples import run_square_plastps

result = run_square_plastps(plot=True)
print("Plastic strain:", result["E"].max())
```

### Plane strain plasticity

```python
from femlabpy.examples import run_hole_plastpe

result = run_hole_plastpe(plot=True)
```

### Custom FEM assembly (low-level)

```python
import femlabpy as fp
import numpy as np

# Initialize arrays
nn = 4  # number of nodes
dof = 2  # degrees of freedom per node
K, p, q = fp.init(nn, dof)

# Element connectivity and coordinates
T = np.array([[1, 2, 3, 4, 1]])  # Q4 element
X = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
G = np.array([[210000, 0.3, 1, 1]])

# Assemble element stiffness
K = fp.kq4e(K, T, X, G)

# Apply boundary conditions
C = np.array([[1, 1, 0], [1, 2, 0], [4, 1, 0]])
P = np.array([[2, 1, 1000], [3, 1, 1000]])

p = fp.setload(p, P)
K, p, _ = fp.setbc(K, p, C, dof)

# Solve
u = np.linalg.solve(K, p)
print("Displacements:", u.ravel())
```

## API Reference

<details>
<summary><strong>Data Loaders</strong></summary>

| Function | Description |
| --- | --- |
| `canti()` | Return cantilever beam benchmark data. Returns dict with T (connectivity), X (coordinates), G (material), C (constraints), P (loads). |
| `flow()` | Return potential flow benchmark data for Q4 and T3 meshes. |
| `bar01()` | Return 2-bar truss benchmark for nonlinear analysis. |
| `bar02()` | Return 3-bar truss benchmark for snap-through analysis. |
| `bar03()` | Return 12-bar dome truss benchmark. |
| `square()` | Return square plate benchmark for plasticity analysis. |
| `hole()` | Return plate with hole benchmark for plasticity analysis. |

</details>

<details>
<summary><strong>High-Level Solvers</strong></summary>

| Function | Description |
| --- | --- |
| `elastic(T, X, G, C, P, dof=2, etype=None)` | Solve a linear elastic T3 or Q4 problem. The wrapper auto-detects the element family from `T` unless `etype` is given explicitly. |
| `flowq4(plot=False)` | Solve potential/thermal problem on Q4 mesh. Returns nodal temperatures and fluxes. |
| `flowt3(plot=False)` | Solve potential/thermal problem on T3 mesh. Returns nodal temperatures and fluxes. |
| `nlbar(T, X, G, C, P, *, no_loadsteps, i_max, i_d, plotdof, tol=1e-8, ...)` | Solve a geometrically nonlinear truss with the legacy orthogonal-residual controls. |
| `plastps(T, X, G, C, P, *, no_loadsteps, i_max, i_d, plotdof, ...)` | Solve the plane-stress elastoplastic legacy driver. |
| `plastpe(T, X, G, C, P, *, no_loadsteps, i_max, i_d, plotdof, ...)` | Solve the plane-strain elastoplastic legacy driver. |

</details>

<details>
<summary><strong>Element Stiffness (ke = element, k = assembly)</strong></summary>

| Function | Description |
| --- | --- |
| `ket3e(Xe, Ge)` | Compute 6x6 stiffness matrix for CST triangle (T3). Xe: 3x2 coordinates, Ge: material [E, nu, t, flag]. |
| `kt3e(K, T, X, G)` | Assemble all T3 element stiffness matrices into global K. |
| `keq4e(Xe, Ge)` | Compute 8x8 stiffness matrix for bilinear quad (Q4) using 2x2 Gauss integration. |
| `kq4e(K, T, X, G)` | Assemble all Q4 element stiffness matrices into global K. |
| `kebar(Xe, A, E, u_e)` | Compute 4x4 tangent stiffness for geometrically nonlinear bar. Includes geometric stiffness. |
| `kbar(K, T, X, G, u)` | Assemble all bar tangent stiffness matrices into global K. |
| `keT4e(Xe, Ge)` | Compute 12x12 stiffness matrix for 4-node tetrahedron (T4). |
| `kT4e(K, T, X, G)` | Assemble all T4 element stiffness matrices into global K. |
| `keh8e(Xe, Ge)` | Compute 24x24 stiffness matrix for 8-node hexahedron (H8) using 2x2x2 Gauss integration. |
| `kh8e(K, T, X, G)` | Assemble all H8 element stiffness matrices into global K. |
| `ket3p(Xe, Ge)` | Compute 3x3 conductivity matrix for T3 potential element. |
| `kt3p(K, T, X, G)` | Assemble all T3 potential element matrices into global K. |
| `keq4p(Xe, Ge)` | Compute 4x4 conductivity matrix for Q4 potential element. |
| `kq4p(K, T, X, G)` | Assemble all Q4 potential element matrices into global K. |
| `keq4eps(Xe, Ge, u_e, H)` | Compute consistent tangent stiffness for plane stress plastic Q4. |
| `kq4eps(K, T, X, G, u, H)` | Assemble plane stress plastic Q4 tangent matrices. |
| `keq4epe(Xe, Ge, u_e, H)` | Compute consistent tangent stiffness for plane strain plastic Q4. |
| `kq4epe(K, T, X, G, u, H)` | Assemble plane strain plastic Q4 tangent matrices. |

</details>

<details>
<summary><strong>Element Response (qe = element, q = assembly)</strong></summary>

| Function | Description |
| --- | --- |
| `qet3e(Xe, Ge, u_e)` | Compute stress and strain for single T3 element. Returns (stress, strain). |
| `qt3e(q, T, X, G, u)` | Compute T3 stresses for all elements. Returns internal forces and stress/strain arrays. |
| `qeq4e(Xe, Ge, u_e)` | Compute stress and strain at 4 Gauss points for single Q4. Returns 4x3 arrays. |
| `qq4e(q, T, X, G, u)` | Compute Q4 stresses for all elements and assemble internal forces. |
| `qebar(Xe, A, E, u_e)` | Compute internal force for single nonlinear bar. Returns axial force and strain. |
| `qbar(q, T, X, G, u)` | Compute bar internal forces and assemble into global vector. |
| `qeT4e(Xe, Ge, u_e)` | Compute stress and strain for single T4 element. |
| `qT4e(q, T, X, G, u)` | Compute T4 stresses for all elements and assemble internal forces. |
| `qeh8e(Xe, Ge, u_e)` | Compute stress and strain at 8 Gauss points for single H8. |
| `qh8e(q, T, X, G, u)` | Compute H8 stresses for all elements and assemble internal forces. |
| `qeq4eps(Xe, Ge, u_e, H)` | Update plane stress plastic Q4 response. Returns stress, strain, updated history. |
| `qq4eps(q, T, X, G, u, H)` | Compute plane stress Q4 internal forces with plasticity. |
| `qeq4epe(Xe, Ge, u_e, H)` | Update plane strain plastic Q4 response. Returns stress, strain, updated history. |
| `qq4epe(q, T, X, G, u, H)` | Compute plane strain Q4 internal forces with plasticity. |

</details>

<details>
<summary><strong>Assembly and Boundary Conditions</strong></summary>

| Function | Description |
| --- | --- |
| `init(nn, dof, *, dynamic=False, use_sparse=None)` | Initialize the global FEM arrays. Returns `(K, p, q)` by default and `(K, M, p, q)` when `dynamic=True`. |
| `assmk(K, ke, nodes, dof)` | Assemble single element stiffness ke into global K at specified nodes. |
| `assmq(q, qe, nodes, dof)` | Assemble single element force qe into global internal force vector q. |
| `setload(p, P)` | Set nodal loads from the FemLab load table. Replaces existing values at the referenced global DOFs. |
| `addload(p, P)` | Add nodal loads from the FemLab load table. Accumulates with existing values. |
| `setbc(K, p, C, dof)` | Apply Dirichlet BCs using direct elimination. Zeros rows/columns, sets diagonal. |
| `solve_lag(K, p, C, dof)` | Solve with Lagrange multipliers for Dirichlet constraints. |
| `solve_lag_general(K, p, G, Q)` | Solve with general linear constraints Gu = Q. |
| `reaction(q, C, dof, comp=None)` | Extract support reactions from the assembled internal-force vector at constrained DOFs. |
| `rnorm(r, C, dof)` | Compute residual norm excluding constrained DOFs. |

</details>

<details>
<summary><strong>Material Models</strong></summary>

| Function | Description |
| --- | --- |
| `devstress(S)` | Compute deviatoric stress and mean stress. S: [sxx, syy, sxy] or [sxx, syy, szz, sxy, syz, sxz]. |
| `eqstress(S)` | Compute von Mises equivalent stress from stress vector. |
| `yieldvm(S, mat, ep, dL)` | Evaluate von Mises yield function f = seq - (sy + H*ep). mat: [E, nu, sy, H]. |
| `dyieldvm(S, mat, ep, dL)` | Derivative of yield function with respect to plastic multiplier dL. |
| `stressvm(S_trial, mat, ep)` | Perform von Mises radial return mapping. Returns (S_updated, delta_ep). |
| `stressdp(S_trial, mat, ep)` | Perform Drucker-Prager return mapping with Newton iterations. |

</details>

<details>
<summary><strong>Mesh I/O</strong></summary>

| Function | Description |
| --- | --- |
| `load_gmsh(filename)` | Load Gmsh mesh file (.msh). Returns GmshMesh with positions, triangles, quads, etc. Legacy 2.x format. |
| `load_gmsh2(filename)` | Load Gmsh mesh with flexible format detection. Supports both 2.x and 4.x (requires [mesh] extra). |

GmshMesh attributes:
- `positions`: Nx3 node coordinates
- `triangles`: Mx4 triangle connectivity [elem_id, n1, n2, n3]
- `quads`: Mx5 quad connectivity
- `lines`: Mx3 line connectivity
- `bounds_min`, `bounds_max`: Bounding box

</details>

<details>
<summary><strong>Plotting</strong></summary>

| Function | Description |
| --- | --- |
| `plotelem(T, X, line_style="k-", nonum=False, noelem=False)` | Plot the undeformed mesh and optionally annotate node or element ids. |
| `plotforces(T, X, P)` | Plot load arrows on the mesh. |
| `plotbc(T, X, C)` | Plot boundary-condition markers. |
| `plotu(T, X, u, dof=None, component=0)` | Plot a scalar nodal field or a contour extracted from a flattened displacement vector. `component=0` plots magnitude, `1` plots x, `2` plots y. |
| `plotq4(T, X, S, scomp)` | Plot a one-based Q4 stress or strain component reconstructed from Gauss-point data. |
| `plott3(T, X, S, scomp)` | Plot a one-based T3 stress or strain component. |

</details>

<details>
<summary><strong>Examples Module</strong></summary>

| Function | Description |
| --- | --- |
| `run_cantilever(plot=False)` | Solve cantilever beam. Returns u, S, E, R. |
| `run_flow_q4(plot=False)` | Solve Q4 potential flow. Returns u (temperatures). |
| `run_flow_t3(plot=False)` | Solve T3 potential flow. Returns u (temperatures). |
| `run_bar01_nlbar(plot=False)` | Solve 2-bar snap-through. Returns load-displacement path. |
| `run_bar02_nlbar(plot=False)` | Solve 3-bar snap-through. Returns load-displacement path. |
| `run_bar03_nlbar(plot=False)` | Solve 12-bar dome. Returns load-displacement path. |
| `run_square_plastps(plot=False)` | Solve square plate plane stress plasticity. |
| `run_square_plastpe(plot=False)` | Solve square plate plane strain plasticity. |
| `run_hole_plastps(plot=False)` | Solve plate with hole plane stress plasticity. |
| `run_hole_plastpe(plot=False)` | Solve plate with hole plane strain plasticity. |
| `run_ex_lag_mult(plot=False)` | Solve 3-bar truss with displacement constraint. |
| `run_gmsh_triangle(plot=False)` | Solve imported Gmsh triangle mesh. |

Data loaders: `cantilever_data()`, `flow_data()`, `bar01_data()`, `bar02_data()`, `bar03_data()`, `square_data(plane_strain=False)`, `hole_data(plane_strain=False)`, `ex_lag_mult_data()`, `gmsh_triangle_data()`

</details>

<details>
<summary><strong>Dynamic Analysis</strong></summary>

| Function | Description |
| --- | --- |
| `solve_newmark(M, C, K, p_func, u0, v0, dt, nsteps, ...)` | Solve time history using implicit Newmark-beta method. Returns `TimeHistory` object. |
| `solve_central_diff(M_lump, C, K, p_func, u0, v0, dt, ...)` | Solve time history using explicit central difference method. Requires lumped mass. |
| `solve_hht(M, C, K, p_func, u0, v0, dt, nsteps, ...)` | Solve time history using HHT-alpha method for high-frequency numerical damping. |
| `solve_modal(K, M, n_modes, C_bc, dof)` | Compute natural frequencies and mode shapes. Returns `ModalResult` object. |
| `seismic_load(M, direction, accel_record, dt_record)` | Create an effective seismic load function $p(t) = -M a_g(t)$ from ground acceleration. |
| `compute_frf(M, C, K, input_dof, output_dof, freq_range)` | Compute Frequency Response Function (FRF) for harmonic excitation. |
| `rayleigh_damping(M, K, f1, f2, zeta1, zeta2)` | Compute Rayleigh damping matrix $C = \alpha M + \beta K$ from two modal frequencies. |
| `plot_time_history(result, dof_index)` | Plot displacement/velocity/acceleration vs time. |
| `plot_modes(X, result, scale)` | Plot deformed mode shapes. |

</details>

<details>
<summary><strong>Periodic Boundaries</strong></summary>

| Function | Description |
| --- | --- |
| `find_periodic_pairs(X, n1_nodes, n2_nodes, tol)` | Find matching node pairs on opposite boundaries for periodic conditions. |
| `periodic_constraints(pairs, dof)` | Generate linear constraint equations (G matrix) for periodic boundary conditions. |
| `homogenize(T, X, G, pairs, dof)` | Compute the 3x3 effective homogenized stiffness tensor (C_eff) for a unit cell. |

</details>

## Dynamic Analysis & Benchmarks

Dynamic analysis examples and benchmarks against OpenSeesPy and CalculiX.

### 1. Static and Modal Analysis (Cantilever Beam)

2D plane-stress cantilever beam (32x4 Q4 elements). Tested for static tip deflection and natural frequencies.

| Metric | Analytical (E-B) | OpenSeesPy | CalculiX (CPS4) | femlabpy (Q4) | femlabpy vs OpenSees |
|---|---|---|---|---|---|
| **Static Deflection** | -9.752e-05 m | -9.544e-05 m | -9.518e-05 m | **-9.544e-05 m** | **0.000%** (Exact) |
| **Frequency: Mode 1** | 26.110 Hz | 26.228 Hz | 26.283 Hz | **26.248 Hz** | **0.07%** |
| **Frequency: Mode 2** | 163.628 Hz | 153.716 Hz | 154.668 Hz | **154.437 Hz** | **0.47%** |
| **Frequency: Mode 3** | 458.163 Hz | 323.810 Hz | 323.892 Hz | **323.881 Hz** | **0.02%** |

*Note: Small differences against 1D Euler-Bernoulli theory are expected for 2D plane-stress elements.*

### 2. Seismic Time History Analysis

2D concrete column (4x32 Q4 elements) under horizontal earthquake (Düzce, BOL090.AT2, PGA = 0.82g). Solved with Newmark-beta average acceleration and 5% Rayleigh damping for 5590 steps (dt = 0.01s).

| Metric | OpenSeesPy (UniformExcitation) | CalculiX (*DLOAD GRAV) | femlabpy (seismic_load) |
|---|---|---|---|
| **Max Roof Disp. (X)** | 5.850 mm | 4.180 mm | **5.848 mm** |
| **Solve Time** | 1.33 seconds | 141.37 seconds | **1.12 seconds** |

**Conclusion:** femlabpy matches OpenSeesPy with a 0.67% RMS difference over 5590 steps. femlabpy solve time is 1.12 seconds.

## Development

```bash
pytest -q                    # run tests
python -m build              # build package
```

## Links

- Source: <https://github.com/adzetto/femlabpy>
- PyPI: <https://pypi.org/project/femlabpy/>
- Issues: <https://github.com/adzetto/femlabpy/issues>

## References

1. Bathe, K.J. (2014). *Finite Element Procedures*. Prentice Hall.
2. Zienkiewicz, Taylor, Zhu (2013). *The Finite Element Method*. Elsevier.
3. Hughes, T.J.R. (2000). *The Finite Element Method*. Dover.
4. Cook et al. (2002). *Concepts and Applications of FEA*. Wiley.
5. de Souza Neto et al. (2008). *Computational Methods for Plasticity*. Wiley.
6. Gmsh: <https://gmsh.info/>
