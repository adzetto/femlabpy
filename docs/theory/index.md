# Theory and Reference Manual

::::{container} fp-home-intro
:::{container} fp-kicker
Implementation-aligned derivations
:::

:::{container} fp-lead
These chapters explain the mathematical structure behind the actual `femlabpy`
implementation. The goal is not to repeat generic FEM lecture notes, but to
connect equations, indexing rules, and constitutive updates to the exact array
operations used in the source tree.
:::
::::

::::{container} fp-card-grid
:::{container} fp-card
## Scope

- Core assembly and boundary-condition logic.
- Bar, triangle, quadrilateral, and solid element formulations.
- Plasticity, modal analysis, damping, dynamic time integration, and periodicity.
- Mesh I/O, plotting, and legacy wrapper behavior.
:::

:::{container} fp-card
## Suggested Order

1. Start with the assembly and boundary chapters.
2. Continue into the element-family chapters for the formulations you use.
3. Read the later chapters for dynamics, periodicity, nonlinear solvers, and I/O.
:::
::::

## Chapter Groups

::::{container} fp-card-grid
:::{container} fp-card
### Core Mechanics

- **01:** Assembly logic and global indexing.
- **02:** Boundary conditions and load treatment.
- **03:** One-dimensional bar formulations.
:::

:::{container} fp-card
### Element Families

- **04:** Three-node triangular elements.
- **05:** Quadrilateral elements, including elastoplastic variants.
- **06:** Tetrahedral and hexahedral solids.
- **07:** Stress invariants and constitutive return mapping.
:::

:::{container} fp-card
### Dynamics, Periodicity, and I/O

- **08:** Dynamics and time integration.
- **09:** Modal reduction and periodic constraints.
- **10:** Mesh import and normalized data structures.
- **11:** Nonlinear legacy solver drivers.
- **12:** Plotting and post-processing.
- **13:** Legacy wrappers and compatibility behavior.
:::
::::

:::{container} fp-section-note
Section titles in this manual are kept as plain text on purpose. Mathematical
symbols and derivations stay in the body text so the sidebar and section
navigation remain clean and readable.
:::

```{toctree}
:maxdepth: 2

01_core_assembly
02_boundary_loads
03_1d_bars
04_2d_triangles
05_2d_quads
06_3d_solids
07_plasticity
08_dynamics
09_modal_periodic
10_io_mesh
11_nonlinear_solvers
12_plotting_postprocessing
13_legacy_wrappers
```
