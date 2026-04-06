# User Guide

<div class="fp-hide-title"></div>

<div class="fp-hero fp-hero--compact">
  <div class="fp-hero-title">User Guide</div>
  <div class="fp-hero-subtitle">
    Practical reading path for model setup, assembly, materials, dynamics, periodic workflows, and extension points.
  </div>
  <div class="fp-hero-actions">
    <a href="ch01_fundamentals.html" class="fp-btn fp-btn-primary">Start Reading Chapter 1 &rarr;</a>
  </div>
</div>

## Chapter Map

<div class="fp-card-grid fp-card-grid--primary">
  <div class="fp-card fp-card--accent-blue">
    <h2><span class="fp-num-badge">1</span> Setup &amp; Core Tables</h2>
    <span class="fp-chip">2 chapters</span>
    <p>Understand the raw <code>X</code>, <code>T</code>, <code>G</code>, <code>C</code>, and <code>P</code> tables. Learn the indexing rules and the basic solve sequence.</p>
    <ul class="fp-chapter-list">
      <li><a href="ch01_fundamentals.html">Ch 1: Fundamentals</a></li>
      <li><a href="ch02_elements.html">Ch 2: Elements</a></li>
    </ul>
  </div>
  <div class="fp-card fp-card--accent-red">
    <h2><span class="fp-num-badge">2</span> Assembly &amp; Materials</h2>
    <span class="fp-chip">2 chapters</span>
    <p>Connect element kernels to global matrices. Define materials, boundary conditions, and loads.</p>
    <ul class="fp-chapter-list">
      <li><a href="ch03_assembly_statics.html">Ch 3: Assembly &amp; Statics</a></li>
      <li><a href="ch04_materials.html">Ch 4: Materials</a></li>
    </ul>
  </div>
  <div class="fp-card fp-card--accent-green">
    <h2><span class="fp-num-badge">3</span> Dynamics &amp; Extensions</h2>
    <span class="fp-chip">4 chapters</span>
    <p>Advance to modal analysis, time integration, periodicity, and building custom element kernels.</p>
    <ul class="fp-chapter-list">
      <li><a href="ch05_dynamics.html">Ch 5: Dynamics</a></li>
      <li><a href="ch06_periodic_io.html">Ch 6: Periodic I/O</a></li>
      <li><a href="ch08_custom_elements.html">Ch 8: Custom Elements</a></li>
      <li><a href="ch09_dynamic_workflows.html">Ch 9: Dynamic Workflows</a></li>
    </ul>
  </div>
</div>

<div class="fp-callout-grid">
  <div class="fp-callout fp-callout--blue">
    <div class="fp-callout-icon">
      <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="10"/><line x1="12" y1="16" x2="12" y2="12"/><line x1="12" y1="8" x2="12.01" y2="8"/></svg>
    </div>
    <div class="fp-callout-body">
      <strong>Tip</strong>
      <p>Keep <code>src/femlabpy/__init__.py</code> open while reading &mdash; it shows which functions are exposed at package level.</p>
    </div>
  </div>
</div>

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