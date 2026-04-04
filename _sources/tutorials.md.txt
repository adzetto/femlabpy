# Tutorials

Learn how to use **femlabpy** for various finite element tasks, from simple statics to complex nonlinear time-history dynamics.

## 1. Linear Static Analysis (Cantilever Beam)

```python
import numpy as np
import femlabpy as fp

# 1. Load Gmsh geometry
mesh = fp.load_gmsh2("beam.msh")
T, X = mesh.quads.astype(int), mesh.positions[:, :2]

# 2. Material (E, nu, plane_stress, thickness, rho)
G = np.array([[210e9, 0.3, 1, 0.1, 7850]])

# 3. Assemble Stiffness
nn, dof = X.shape[0], 2
K, p = fp.init(nn, dof)
K = fp.kq4e(K, T, X, G)

# 4. Apply BCs and Loads, then Solve
# ... (see examples folder for full scripts)
```

## 2. Dynamic Analysis

```python
from femlabpy import solve_newmark
from femlabpy.dynamics import seismic_load

# Assemble Mass
M = fp.init(nn, dof)[0]
M = fp.mq4e(M, T, X, G)

# Create seismic effective load
direction = np.array([1, 0] * nn)
p_eff = seismic_load(M, direction, accel_data, dt)

# Solve using Newmark-beta
history = solve_newmark(M, np.zeros_like(M), K, p_eff, u0, v0, dt, nsteps=1000)
```
