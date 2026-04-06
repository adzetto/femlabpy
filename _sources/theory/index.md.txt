# Theory and Reference Manual

<div class="fp-hero fp-hero--compact">
  <div class="fp-hero-title">Theory &amp; Reference Manual</div>
  <div class="fp-hero-subtitle">
    Implementation-aligned derivations that connect equations, indexing rules, and constitutive updates to the exact array operations used in the source tree.
  </div>
</div>

<div class="fp-card-grid fp-card-grid--primary">
  <div class="fp-card fp-card--accent-blue">
    <h3>Core Mechanics</h3>
    <span class="fp-chip">3 chapters</span>
    <ul class="fp-chapter-list">
      <li><a href="01_core_assembly.html">01 &mdash; Assembly &amp; Global Indexing</a></li>
      <li><a href="02_boundary_loads.html">02 &mdash; Boundary Conditions &amp; Loads</a></li>
      <li><a href="03_1d_bars.html">03 &mdash; 1-D Bar Formulations</a></li>
    </ul>
  </div>
  <div class="fp-card fp-card--accent-red">
    <h3>Element Families</h3>
    <span class="fp-chip">4 chapters</span>
    <ul class="fp-chapter-list">
      <li><a href="04_2d_triangles.html">04 &mdash; Triangular Elements (T3)</a></li>
      <li><a href="05_2d_quads.html">05 &mdash; Quadrilateral Elements (Q4)</a></li>
      <li><a href="06_3d_solids.html">06 &mdash; Tetrahedra &amp; Hexahedra</a></li>
      <li><a href="07_plasticity.html">07 &mdash; Constitutive Return Mapping</a></li>
    </ul>
  </div>
  <div class="fp-card fp-card--accent-green">
    <h3>Dynamics, Periodicity &amp; I/O</h3>
    <span class="fp-chip">6 chapters</span>
    <ul class="fp-chapter-list">
      <li><a href="08_dynamics.html">08 &mdash; Time Integration</a></li>
      <li><a href="09_modal_periodic.html">09 &mdash; Modal &amp; Periodic</a></li>
      <li><a href="10_io_mesh.html">10 &mdash; Mesh Import</a></li>
      <li><a href="11_nonlinear_solvers.html">11 &mdash; Nonlinear Solvers</a></li>
      <li><a href="12_plotting_postprocessing.html">12 &mdash; Plotting</a></li>
      <li><a href="13_legacy_wrappers.html">13 &mdash; Legacy Wrappers</a></li>
    </ul>
  </div>
</div>

<div class="fp-callout-grid">
  <div class="fp-callout fp-callout--blue">
    <div class="fp-callout-icon">
      <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M2 3h6a4 4 0 014 4v14a3 3 0 00-3-3H2z"/><path d="M22 3h-6a4 4 0 00-4 4v14a3 3 0 013-3h7z"/></svg>
    </div>
    <div class="fp-callout-body">
      <strong>Reading order</strong>
      <p>Start with assembly (01-02), then the element families you use, then dynamics and I/O last.</p>
    </div>
  </div>
</div>

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
