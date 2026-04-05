# User Guide

<div class="fp-hide-title"></div>

<div class="fp-hero" style="padding: 3rem 1.5rem; margin-bottom: 3rem; background: linear-gradient(180deg, var(--fp-surface) 0%, transparent 100%);">
  <div class="fp-hero-title" style="font-size: 3rem;">User Guide</div>
  <div class="fp-hero-subtitle">
    Practical reading path for model setup, assembly, materials, dynamics, periodic workflows, and extension points.
  </div>
  <div class="fp-hero-actions">
    <a href="ch01_fundamentals.html" class="fp-btn fp-btn-primary">Start Reading Chapter 1 →</a>
  </div>
</div>

## Chapter Map

<div class="fp-card-grid">
  <div class="fp-card">
    <h2><span style="font-size: 1.5em; vertical-align: middle;">1️⃣</span> Setup & Core Tables</h2>
    <p>Understand the raw <code>X</code>, <code>T</code>, <code>G</code>, <code>C</code>, and <code>P</code> tables. Learn the indexing rules and the basic solve sequence.</p>
    <ul style="margin-top: 1rem; color: var(--fp-accent); list-style: none; padding: 0; font-weight: 500;">
      <li style="margin-bottom: 0.5rem;"><a href="ch01_fundamentals.html">Chapter 1: Fundamentals</a></li>
      <li><a href="ch02_elements.html">Chapter 2: Elements</a></li>
    </ul>
  </div>
  <div class="fp-card">
    <h2><span style="font-size: 1.5em; vertical-align: middle;">2️⃣</span> Assembly & Materials</h2>
    <p>Connect element kernels to global matrices. Define materials, boundary conditions, and loads.</p>
    <ul style="margin-top: 1rem; color: var(--fp-accent); list-style: none; padding: 0; font-weight: 500;">
      <li style="margin-bottom: 0.5rem;"><a href="ch03_assembly_statics.html">Chapter 3: Assembly & Statics</a></li>
      <li><a href="ch04_materials.html">Chapter 4: Materials</a></li>
    </ul>
  </div>
  <div class="fp-card">
    <h2><span style="font-size: 1.5em; vertical-align: middle;">3️⃣</span> Dynamics & Extensions</h2>
    <p>Advance to modal analysis, time integration, periodicity, and building custom element kernels.</p>
    <ul style="margin-top: 1rem; color: var(--fp-accent); list-style: none; padding: 0; font-weight: 500;">
      <li style="margin-bottom: 0.5rem;"><a href="ch05_dynamics.html">Chapter 5: Dynamics</a></li>
      <li style="margin-bottom: 0.5rem;"><a href="ch06_periodic_io.html">Chapter 6: Periodic I/O</a></li>
      <li style="margin-bottom: 0.5rem;"><a href="ch08_custom_elements.html">Chapter 8: Custom Elements</a></li>
      <li><a href="ch09_dynamic_workflows.html">Chapter 9: Dynamic Workflows</a></li>
    </ul>
  </div>
</div>

:::{container} fp-section-note
If you are reading the source code while using this manual, keep
`src/femlabpy/__init__.py` open in another tab. It shows which functions are
exposed at package level and helps you connect the public API to the underlying
implementation modules.
:::

```{toctree}
:maxdepth: 2
:hidden:

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