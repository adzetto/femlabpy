(index)=
# femlabpy

<div class="fp-hide-title"></div>

<div class="fp-hero">
  <div class="fp-hero-title">femlabpy</div>
  <div class="fp-hero-subtitle">
    A modern finite element teaching and research library built around explicit NumPy arrays, compact element kernels, and source-level transparency.
  </div>
  <div class="fp-hero-actions">
    <a href="tutorials.html" class="fp-btn fp-btn-primary">Get Started with Tutorials</a>
    <a href="manual/index.html" class="fp-btn fp-btn-secondary">Read the User Guide</a>
  </div>
</div>

```{toctree}
:maxdepth: 2
:hidden:

Tutorials <tutorials>
Guide <manual/index>
Theory <theory/index>
API <api>
MATLAB Mapping <matlab_python_mapping>
Roadmap & Arch <roadmap>
```

<div class="fp-card-grid">
  <div class="fp-card">
    <h2><span style="font-size: 1.5em; vertical-align: middle;">📚</span> User Guide</h2>
    <p>Practical, end-to-end reading for model setup, assembly, materials, dynamics, periodic workflows, and extension points.</p>
    <a href="manual/index.html">Open the guide →</a>
  </div>
  <div class="fp-card">
    <h2><span style="font-size: 1.5em; vertical-align: middle;">📐</span> Theory and Reference</h2>
    <p>Implementation-aligned derivations for assembly, elements, plasticity, time integration, periodic constraints, and mesh import.</p>
    <a href="theory/index.html">Open the theory manual →</a>
  </div>
  <div class="fp-card">
    <h2><span style="font-size: 1.5em; vertical-align: middle;">⚙️</span> API Reference</h2>
    <p>Module-level documentation that explains ownership, reading order, and the full function and class reference generated from the source code.</p>
    <a href="api.html">Open the API reference →</a>
  </div>
  <div class="fp-card">
    <h2><span style="font-size: 1.5em; vertical-align: middle;">🔄</span> MATLAB Mapping</h2>
    <p>Cross-reference notes for readers moving from the original classroom scripts to the Python port and its public API.</p>
    <a href="matlab_python_mapping.html">Open the mapping notes →</a>
  </div>
  <div class="fp-card">
    <h2><span style="font-size: 1.5em; vertical-align: middle;">🚀</span> Roadmap & Architecture</h2>
    <p>View the implementation plans, development roadmap, and deep architectural comparisons with other solvers like OpenSees.</p>
    <a href="roadmap.html">View the roadmap →</a>
  </div>
</div>

## What This Documentation Covers

<div class="fp-stat-grid">
  <div class="fp-stat">
    <strong>9 manual chapters</strong>
    <p>Practical workflows from array layout to dynamic post-processing.</p>
  </div>
  <div class="fp-stat">
    <strong>13 theory chapters</strong>
    <p>Derivations and implementation notes matched to the actual source tree.</p>
  </div>
  <div class="fp-stat">
    <strong>18 API modules</strong>
    <p>Core workflow, element kernels, dynamics, periodicity, plotting, and mesh I/O.</p>
  </div>
  <div class="fp-stat">
    <strong>Static to transient</strong>
    <p>Linear statics, modal analysis, periodic RVEs, and Newmark-family solvers.</p>
  </div>
</div>

## Reading Order

<div class="fp-card-grid">
  <div class="fp-card">
    <h3>1. Start With Tutorials</h3>
    <p>Read the <a href="tutorials.html">tutorials</a> for a compact orientation before diving into the module pages.</p>
  </div>
  <div class="fp-card">
    <h3>2. Learn The Workflow</h3>
    <p>Use the <a href="manual/index.html">User Guide</a> to understand how `X`, `T`, `G`, `C`, and `P` move through assembly and solution.</p>
  </div>
  <div class="fp-card">
    <h3>3. Match Equations To Code</h3>
    <p>Use the <a href="theory/index.html">Theory manual</a> when you want the mathematical meaning behind the actual vectorized implementation.</p>
  </div>
  <div class="fp-card">
    <h3>4. Look Up Exact Behavior</h3>
    <p>Use the <a href="api.html">API reference</a> when you need precise call signatures, return types, module ownership, and the full reference text.</p>
  </div>
</div>

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