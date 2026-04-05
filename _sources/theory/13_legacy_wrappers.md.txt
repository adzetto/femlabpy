---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Chapter 13: Legacy Wrappers

This chapter explains the compatibility layer in `matlab.py` and `compat.py`.
The goal is not to hide the modern Python API. The goal is to let older FemLab
scripts keep working while the underlying code uses the shared core modules.

## Why The Wrappers Exist

The original FemLab material was organized as single-call drivers such as
`elastic`, `flowq4`, `nlbar`, `plastps`, and `plastpe`. The wrapper layer keeps
that entry point available, but each wrapper delegates to the shared core
routines rather than duplicating the implementation.

That separation matters:

- the wrappers preserve the old classroom interface,
- the core modules stay reusable and testable,
- the plotting and recovery logic remains centralized.

## Shared Helpers

### `_legacy_bundle`

`_legacy_bundle(T, X, G, C, P, dof=None)` converts the legacy input deck into a
typed dictionary with normalized arrays. It resolves the degree-of-freedom count
from `X` when needed and creates an empty load table when `P` is omitted.

The wrappers use this helper so every legacy driver sees the same normalized
shape and dtype rules.

### `_set_axis`

`_set_axis(ax, limits)` applies MATLAB-style axis vectors to Matplotlib axes.
It supports both 2D and 3D limit arrays and is used only for compatibility with
the old plotting scripts.

### `compat.setpath`

`compat.setpath()` returns the package and examples directories and can append
the examples directory to `sys.path`.

That mirrors the original MATLAB `setpath.m` helper closely enough for old
scripts that expect the examples directory to be discoverable without manual
path setup.

## Benchmark Data Wrappers

The simplest wrappers just return packaged benchmark decks:

- `canti()`
- `flow()`
- `bar01()`
- `bar02()`
- `bar03()`
- `square()`
- `hole()`

These functions are data accessors, not solvers. They give you the original
coordinates, topology, materials, loads, and constraints in the exact format
expected by the legacy drivers.

`square()` and `hole()` both accept `plane_strain` so the caller can choose the
plane-stress or plane-strain benchmark variant.

## Solver Wrappers

### `elastic`

`elastic(T, X, G, C, P, ...)` solves the linear Q4 elasticity workflow. It uses
`init`, `kq4e`, `setload`, `setbc`, `qq4e`, and `reaction` under the hood.

If `plot=True`, it also returns figures for the mesh, loads, constraints, and
stress field. That makes it a faithful compatibility wrapper for the original
single-call MATLAB driver.

### `flowq4` and `flowt3`

These wrappers solve the scalar potential benchmark using either Q4 or T3
elements. Internally they both call `_flow_driver`, which chooses the correct
element kernels:

- `kq4p` and `qq4p` for quadrilaterals,
- `kt3p` and `qt3p` for triangles.

The return value contains the solved field, the recovered gradients, and an
optional figure that matches the legacy classroom output.

### `nlbar`

`nlbar(...)` wraps `solve_nlbar`. It keeps the original load-step controls and
optionally builds the load-displacement and deformed-geometry plots.

### `plastps` and `plastpe`

These wrappers call `_plast_driver`, which in turn uses `solve_plastic`.

- `plastps` configures plane stress.
- `plastpe` configures plane strain.

Both wrappers preserve the old MATLAB-style plot controls and can return the
equivalent plastic strain field in the deformed configuration.

## What The Wrappers Return

The wrappers return dictionaries rather than bare arrays. That is a deliberate
departure from the original MATLAB style because it keeps the response data and
the optional figures together in one object.

The most common keys are:

- `u`
- `q`
- `S`
- `E`
- `R`
- `data`
- `figures`

That structure makes it easier to inspect results programmatically and keeps the
legacy wrappers usable in notebooks and scripts.

## Reading Order

If you need to modify the compatibility layer, read these files in this order:

1. `compat.py`
2. `matlab.py`
3. `core.py`
4. `boundary.py`
5. `loads.py`
6. `elements/`
7. `solvers.py`
8. `plotting.py`
9. `postprocess.py`

That path matches the actual dependency chain used by the wrappers.

