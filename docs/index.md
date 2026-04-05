# femlabpy

::::{container} fp-home-intro
:::{container} fp-kicker
Finite element workflows, theory, and source-level reference
:::

:::{container} fp-lead
`femlabpy` is a finite element teaching and research library built around
explicit NumPy arrays, compact element kernels, and workflows that stay close
to the original FemLab MATLAB and Scilab material while using modern Python
tooling.
:::
::::

```{toctree}
:maxdepth: 2
:hidden:

Tutorials <tutorials>
Guide <manual/index>
Theory <theory/index>
API <api>
MATLAB Mapping <matlab_python_mapping>
```

::::{container} fp-card-grid
:::{container} fp-card
## User Guide

Practical, end-to-end reading for model setup, assembly, materials, dynamics,
periodic workflows, and extension points.

[Open the guide](manual/index.md)
:::

:::{container} fp-card
## Theory and Reference

Implementation-aligned derivations for assembly, elements, plasticity,
time integration, periodic constraints, and mesh import.

[Open the theory manual](theory/index.md)
:::

:::{container} fp-card
## API Reference

Module-level documentation that explains ownership, reading order, and the full
function and class reference generated from the source code.

[Open the API reference](api.rst)
:::

:::{container} fp-card
## MATLAB Mapping

Cross-reference notes for readers moving from the original classroom scripts to
the Python port and its public API.

[Open the mapping notes](matlab_python_mapping.md)
:::
::::

## What This Documentation Covers

::::{container} fp-stat-grid
:::{container} fp-stat
**9 manual chapters**

Practical workflows from array layout to dynamic post-processing.
:::

:::{container} fp-stat
**13 theory chapters**

Derivations and implementation notes matched to the actual source tree.
:::

:::{container} fp-stat
**18 API modules**

Core workflow, element kernels, dynamics, periodicity, plotting, and mesh I/O.
:::

:::{container} fp-stat
**Static to transient**

Linear statics, modal analysis, periodic RVEs, and Newmark-family solvers.
:::
::::

## Reading Order

::::{container} fp-card-grid
:::{container} fp-card
### 1. Start With Tutorials

Read [tutorials](tutorials.md) for a compact orientation before diving into the
module pages.
:::

:::{container} fp-card
### 2. Learn The Workflow

Use the [User Guide](manual/index.md) to understand how `X`, `T`, `G`, `C`,
and `P` move through assembly and solution.
:::

:::{container} fp-card
### 3. Match Equations To Code

Use the [Theory manual](theory/index.md) when you want the mathematical meaning
behind the actual vectorized implementation.
:::

:::{container} fp-card
### 4. Look Up Exact Behavior

Use the [API reference](api.rst) when you need precise call signatures, return
types, module ownership, and the full reference text.
:::
::::

:::{important}
This documentation is intentionally split into layers.

- The guide explains how to use the library in a realistic workflow.
- The theory pages explain why the kernels and solvers are written that way.
- The API reference explains exactly which module owns which part of the code.
:::

:::{note}
Many public tables in `femlabpy` use one-based node or property numbering to
preserve FemLab compatibility, while the underlying NumPy operations remain
zero-based. That convention is repeated throughout the guide and theory pages
because it is one of the main things to understand before reading the kernels.
:::
