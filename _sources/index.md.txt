# femlabpy

Welcome to the **femlabpy** documentation!

femlabpy is a powerful, pure-Python finite element library designed for both educational clarity and production-grade performance. It brings the intuitive API of MATLAB/Scilab FemLab into the modern Python scientific ecosystem.

```{toctree}
:maxdepth: 2
:hidden:

tutorials
api
```

## Quick Installation

```bash
pip install femlabpy
```

## Features

- **Elements:** Truss, Plane Stress/Strain (T3, Q4), 3D Solids (T4, H8).
- **Material Models:** Elastic, von Mises Plasticity, Drucker-Prager.
- **Dynamic Analysis:** Modal extraction, Newmark-beta, HHT-alpha time-history integration.
- **Advanced Tools:** Periodic boundary conditions for homogenization.
- **I/O:** Full `gmsh` support for loading industrial-scale meshes.

Check out the [Tutorials](tutorials.md) or the detailed [API Reference](api.rst) for more!
