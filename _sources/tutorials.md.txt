# Tutorials

This page guides you through using femlabpy. You will learn how to build a mesh, apply loads, and solve finite element problems. We use short, practical code blocks.

## 1. Installation

First, install femlabpy and its mesh tools.

```bash
pip install "femlabpy[mesh]"
```

We recommend installing `gmsh` to create and load geometry.

---

## 2. Basic Static Analysis

In this example, we model a cantilever beam fixed on the left and loaded on the right.

### Step 1: Create a Mesh

You can create a mesh using the Gmsh Python API. We will make a 4x0.5 rectangle.

```python
import gmsh
import numpy as np

gmsh.initialize()
gmsh.model.add("beam")

# Create a 4.0 by 0.5 rectangle
gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 4.0, 0.5)
gmsh.model.occ.synchronize()

# Mesh it with structured quads
gmsh.model.mesh.setRecombine(2, 1)  # 2D surface 1
gmsh.model.mesh.generate(2)
gmsh.write("beam.msh")
gmsh.finalize()
```

### Step 2: Load the Mesh

Use femlabpy to read the `.msh` file.

```python
import femlabpy as fp

mesh = fp.load_gmsh2("beam.msh")

# Extract quad elements and node coordinates
T = mesh.quads.astype(int)
X = mesh.positions[:, :2]

nn = X.shape[0]  # Number of nodes
dof = 2          # Degrees of freedom per node
```

### Step 3: Properties and Assembly

Set the material table `G`. The order for 2D is: Young's Modulus, Poisson's Ratio, Plane Stress (1) or Strain (2), Thickness, Density.

```python
# E=210GPa, nu=0.3, Plane Stress, t=0.1m, rho=7850kg/m3
G = np.array([[210e9, 0.3, 1.0, 0.1, 7850.0]])

# Initialize empty stiffness matrix K and load vector p
K, p = fp.init(nn, dof)

# Assemble global stiffness matrix
K = fp.kq4e(K, T, X, G)
```

### Step 4: Boundary Conditions and Solving

Fix the left edge (x=0) and apply a downward force at the top right node.

```python
# Find nodes on the left edge
tol = 1e-6
left_nodes = np.where(X[:, 0] < tol)[0] + 1  # +1 because femlabpy uses 1-based node IDs

# Create constraint array C: [node_id, dof, value]
C = []
for n in left_nodes:
    C.append([n, 1, 0.0]) # fix X
    C.append([n, 2, 0.0]) # fix Y
C = np.array(C)

# Find top right node
right_nodes = np.where((X[:, 0] > 4.0 - tol) & (X[:, 1] > 0.5 - tol))[0] + 1
P = np.array([[right_nodes[0], 2, -1000.0]]) # -1000 N in Y direction

# Apply loads and constraints
p = fp.setload(p, P, dof)
K_bc, p_bc, _ = fp.setbc(K, p, C, dof)

# Solve the linear system
u = np.linalg.solve(K_bc, p_bc)

# The result 'u' is a 1D array of displacements
print(f"Max displacement: {np.min(u):.6e} m")
```

---

## 3. Dynamic Analysis

You can also solve time-dependent problems. We use the same stiffness matrix `K` from above.

### Modal Analysis

Find the natural frequencies and mode shapes.

```python
# Assemble mass matrix
M = np.zeros_like(K)
M = fp.mq4e(M, T, X, G)

# Solve for the first 3 modes
modal_result = fp.solve_modal(K, M, n_modes=3, C_bc=C, dof=2)
print("Frequencies (Hz):", modal_result.freq_hz)
```

### Time-History Analysis (Earthquake)

Apply a ground acceleration to the base.

```python
from femlabpy.dynamics import seismic_load, solve_newmark
from femlabpy.damping import rayleigh_damping

# Create an effective seismic load for horizontal (X) shaking
# inf_vec has 1.0 for X DOFs and 0.0 for Y DOFs
inf_vec = np.zeros(nn * dof)
inf_vec[0::2] = 1.0

# Assume accel_data is a numpy array of ground acceleration in m/s^2
dt = 0.01
p_eff = seismic_load(M, inf_vec, accel_data, dt)

# Add 5% Rayleigh damping based on the first two modes
from femlabpy.damping import rayleigh_coefficients
a, b = rayleigh_coefficients(modal_result.freq_hz[0], modal_result.freq_hz[1], 0.05, 0.05)
C_damp = rayleigh_damping(M, K, a, b)

# Initial conditions
u0 = np.zeros((nn * dof, 1))
v0 = np.zeros((nn * dof, 1))

# Solve using implicit Newmark-beta
nsteps = len(accel_data)
result = solve_newmark(M, C_damp, K, p_eff, u0, v0, dt, nsteps, C_bc=C)

# result.u contains the displacement history
print("Max dynamic roof displacement:", np.max(np.abs(result.u)))
```

---

## 4. Periodic Boundaries (Homogenization)

If you model a repeating unit cell (like a composite material), you must tie opposite boundaries together so they deform equally.

```python
from femlabpy.periodic import find_periodic_pairs, solve_periodic, homogenize

# 1. Find nodes on the left/right and bottom/top boundaries
left_nodes = np.where(X[:, 0] < tol)[0] + 1
right_nodes = np.where(X[:, 0] > 1.0 - tol)[0] + 1
bottom_nodes = np.where(X[:, 1] < tol)[0] + 1
top_nodes = np.where(X[:, 1] > 1.0 - tol)[0] + 1

# 2. Pair them up
lr_pairs = find_periodic_pairs(X, left_nodes, right_nodes)
bt_pairs = find_periodic_pairs(X, bottom_nodes, top_nodes)
all_pairs = lr_pairs + bt_pairs

# 3. Compute effective stiffness (C_eff) of the unit cell
# This applies 3 macro-strain states and computes the average stress.
C_eff = homogenize(T, X, G, all_pairs, dof=2)

print("Effective Homogenized Modulus Matrix:")
print(C_eff)
```

These simple tools allow you to build scripts for research, testing, or production.
