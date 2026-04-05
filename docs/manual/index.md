# User Guide

::::{container} fp-home-intro
:::{container} fp-kicker
Practical reading path
:::

:::{container} fp-lead
This guide is the hands-on entry point for `femlabpy`. It focuses on how to
prepare arrays, call the library in the right order, and understand what each
stage of the workflow is doing.
:::
::::

::::{container} fp-card-grid
:::{container} fp-card
## Who This Section Is For

- Users coming from MATLAB, Scilab, NumPy, or a course FEM codebase.
- Readers who want implementation-aware tutorials rather than black-box examples.
- Anyone trying to connect the raw `X`, `T`, `G`, `C`, and `P` tables to the
  global matrices and solver outputs.
:::

:::{container} fp-card
## Recommended Path

1. Start with Chapter 1 for the array layout and the global workflow.
2. Read Chapter 2 before touching element kernels directly.
3. Read Chapters 3 and 4 for statics, constraints, and material conventions.
4. Read Chapters 5, 6, and 9 for dynamics, periodicity, and Gmsh-driven work.
5. Read Chapters 7 and 8 when extending the library or adapting examples.
:::
::::

## Chapter Map

::::{container} fp-card-grid
:::{container} fp-card
### Setup and Core Tables

- **Chapter 1:** Core data tables, indexing rules, and the basic solve sequence.
- **Chapter 2:** Element-family overview for bars, triangles, quads, and solids.
:::

:::{container} fp-card
### Assembly and Materials

- **Chapter 3:** Global assembly, loads, supports, and constrained solves.
- **Chapter 4:** Material-table structure and constitutive helper routines.
:::

:::{container} fp-card
### Dynamics and Extensions

- **Chapter 5:** Modal analysis, damping, load callables, and time integration.
- **Chapter 6:** Periodic workflows, homogenization, and Gmsh import.
- **Chapter 7:** Packaged examples and bundled driver organization.
- **Chapter 8:** Custom element development and integration.
- **Chapter 9:** Dynamic response, FRF, seismic loading, and history data.
:::
::::

:::{container} fp-section-note
If you are reading the source code while using this manual, keep
`src/femlabpy/__init__.py` open in another tab. It shows which functions are
exposed at package level and helps you connect the public API to the underlying
implementation modules.
:::

```{toctree}
:maxdepth: 2

ch01_fundamentals
ch02_elements
ch03_assembly_statics
ch04_materials
ch05_dynamics
ch06_periodic_io
ch07_packaged_examples
ch08_custom_elements
ch09_dynamic_workflows
```
