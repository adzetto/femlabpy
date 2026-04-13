(index)=
# femlabpy

<div class="fp-hide-title"></div>

<div class="fp-hero">
  <div class="fp-hero-meta">
    <span class="fp-meta-badge">v0.6.2</span>
    <span class="fp-meta-badge">MIT License</span>
    <span class="fp-meta-badge">Python &ge; 3.9</span>
  </div>
  <div class="fp-hero-subtitle">
    Python finite element toolkit.
  </div>

  <div class="fp-hero-commands">
    <div class="fp-cmd-row">
      <code class="fp-cmd-text" id="cmd-pip">pip install femlabpy</code>
      <button class="fp-cmd-copy" data-target="cmd-pip" aria-label="Copy">
        <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 01-2-2V4a2 2 0 012-2h9a2 2 0 012 2v1"/></svg>
      </button>
    </div>
    <div class="fp-cmd-row">
      <code class="fp-cmd-text" id="cmd-pydoc">python -m pydoc femlabpy</code>
      <button class="fp-cmd-copy" data-target="cmd-pydoc" aria-label="Copy">
        <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 01-2-2V4a2 2 0 012-2h9a2 2 0 012 2v1"/></svg>
      </button>
    </div>
  </div>

  <div class="fp-hero-interactive" id="fp-hero-interactive">
    <div class="fp-naming-anim-container" id="fp-naming-anim">
      <!-- JS will populate the naming animation here -->
    </div>
    <div class="fp-arg-panel" id="fp-arg-panel">
      <!-- JS will populate static arg callouts here -->
    </div>
  </div>

  <div class="fp-hero-actions">
    <a href="theory/index.html" class="fp-btn fp-btn-primary">Explore the Theory</a>
    <a href="manual/index.html" class="fp-btn fp-btn-secondary">Read the User Guide</a>
    <a href="releases.html" class="fp-btn fp-btn-ghost">Browse Releases</a>
    <a href="https://github.com/adzetto/femlabpy" class="fp-btn fp-btn-ghost">View on GitHub</a>
  </div>
</div>

```{toctree}
:maxdepth: 2
:hidden:

Tutorials <tutorials>
Guide <manual/index>
Theory <theory/index>
API <api>
Releases <releases>
MATLAB Mapping <matlab_python_mapping>
Roadmap & Arch <roadmap>
```

<div class="fp-card-grid fp-card-grid--primary">
  <div class="fp-card fp-card--accent-blue">
    <div class="fp-icon"><svg xmlns="http://www.w3.org/2000/svg" height="28" viewBox="0 -960 960 960" width="28"><path d="M560-564v-68q33-14 67.5-21t72.5-7q26 0 51 4t49 10v64q-24-9-48.5-13.5T700-644q-38 0-73 9.5T560-604Zm0 220v-68q33-14 67.5-21t72.5-7q26 0 51 4t49 10v64q-24-9-48.5-13.5T700-424q-38 0-73 9.5T560-384Zm0-110v-68q33-14 67.5-21t72.5-7q26 0 51 4t49 10v64q-24-9-48.5-13.5T700-534q-38 0-73 9.5T560-494ZM260-320q47 0 91.5 10.5T440-278v-394q-41-24-87-36t-93-12q-36 0-71.5 7T120-692v396q35-12 69.5-18t70.5-6Zm260 42q44-21 88.5-31.5T700-320q36 0 70.5 6t69.5 18v-396q-33-14-68.5-21t-71.5-7q-47 0-93 12t-87 36v394Zm-40 118q-48-38-104-59t-116-21q-42 0-82.5 11T100-198q-21 11-40.5-1T40-234v-482q0-11 5.5-21T62-752q46-24 96-36t102-12q58 0 113.5 15T480-740q51-30 106.5-45T700-800q52 0 102 12t96 36q11 5 16.5 15t5.5 21v482q0 23-19.5 35t-40.5 1q-37-20-77.5-31T700-240q-60 0-116 21t-104 59ZM280-494Z"/></svg></div>
    <h2>User Guide</h2>
    <span class="fp-chip">9 chapters</span>
    <p>Practical, end-to-end reading for model setup, assembly, materials, dynamics, periodic workflows, and extension points.</p>
    <a href="manual/index.html">Open the guide &rarr;</a>
  </div>
  <div class="fp-card fp-card--accent-red">
    <div class="fp-icon"><svg xmlns="http://www.w3.org/2000/svg" height="28" viewBox="0 -960 960 960" width="28"><path d="M320-240 80-480l240-240 57 57-184 184 183 183-56 56Zm320 0-57-57 184-184-183-183 56-56 240 240-240 240Z"/></svg></div>
    <h2>Theory &amp; Reference</h2>
    <span class="fp-chip">13 chapters</span>
    <p>Implementation-aligned derivations for assembly, elements, plasticity, time integration, periodic constraints, and mesh import.</p>
    <a href="theory/index.html">Open the theory manual &rarr;</a>
  </div>
  <div class="fp-card fp-card--accent-green">
    <div class="fp-icon"><svg xmlns="http://www.w3.org/2000/svg" height="28" viewBox="0 -960 960 960" width="28"><path d="M320-240h320v-80H320v80Zm0-160h320v-80H320v80ZM240-80q-33 0-56.5-23.5T160-160v-640q0-33 23.5-56.5T240-880h320l240 240v480q0 33-23.5 56.5T720-80H240Zm280-520v-200H240v640h480v-440H520ZM240-800v200-200 640-640Z"/></svg></div>
    <h2>API Reference</h2>
    <span class="fp-chip">18 modules</span>
    <p>Module-level documentation with ownership, reading order, and the full function and class reference from source.</p>
    <a href="api.html">Open the API reference &rarr;</a>
  </div>
</div>

<div class="fp-card-grid fp-card-grid--secondary">
  <div class="fp-card">
    <div class="fp-icon"><svg xmlns="http://www.w3.org/2000/svg" height="28" viewBox="0 -960 960 960" width="28"><path d="M280-80q-33 0-56.5-23.5T200-160v-640h80v640h480v-280h80v280q0 33-23.5 56.5T760-80H280Zm200-240L320-480l56-58 64 64v-326h80v326l64-64 56 58-160 160Zm-200 0v160-160Zm480-200v-280H480v-80h280q33 0 56.5 23.5T840-800v280h-80Z"/></svg></div>
    <h2>Release History</h2>
    <p>Read the tagged release timeline, generated directly from version history and commit ranges, with formulas and raw commit subjects for each version.</p>
    <a href="releases.html">Open the release notes &rarr;</a>
  </div>
  <div class="fp-card">
    <div class="fp-icon"><svg xmlns="http://www.w3.org/2000/svg" height="28" viewBox="0 -960 960 960" width="28"><path d="M160-160v-80h110l-16-14q-52-46-73-105t-21-119q0-111 66.5-197.5T400-790v84q-72 26-116 88.5T240-478q0 45 17 87.5t53 78.5l10 12v-100h80v240H160Zm400-10v-84q72-26 116-88.5T720-482q0-45-17-87.5T650-648l-10-12v100h-80v-240h240v80H690l16 14q49 49 71.5 106.5T800-482q0 111-66.5 197.5T560-170Z"/></svg></div>
    <h2>MATLAB Mapping</h2>
    <p>Cross-reference notes for readers moving from the original classroom MATLAB scripts to the Python port.</p>
    <a href="matlab_python_mapping.html">Open the mapping notes &rarr;</a>
  </div>
  <div class="fp-card">
    <div class="fp-icon"><svg xmlns="http://www.w3.org/2000/svg" height="28" viewBox="0 -960 960 960" width="28"><path d="m136-240-56-56 296-298 160 160 208-206H640v-80h240v240h-80v-104L536-320 376-480 136-240Z"/></svg></div>
    <h2>Roadmap &amp; Architecture</h2>
    <p>Implementation plans, development roadmap, and architectural comparisons with OpenSees and CalculiX.</p>
    <a href="roadmap.html">View the roadmap &rarr;</a>
  </div>
</div>

## Reading Order

<div class="fp-steps">
  <div class="fp-step">
    <span class="fp-step-num">1</span>
    <div class="fp-step-content">
      <h3>Start With Tutorials</h3>
      <p>Read the <a href="tutorials.html">tutorials</a> for a compact orientation before diving into the module pages.</p>
    </div>
  </div>
  <div class="fp-step">
    <span class="fp-step-num">2</span>
    <div class="fp-step-content">
      <h3>Learn The Workflow</h3>
      <p>Use the <a href="manual/index.html">User Guide</a> to understand how <strong>X</strong>, <strong>T</strong>, <strong>G</strong>, <strong>C</strong>, and <strong>P</strong> move through assembly and solution.</p>
    </div>
  </div>
  <div class="fp-step">
    <span class="fp-step-num">3</span>
    <div class="fp-step-content">
      <h3>Match Equations To Code</h3>
      <p>Use the <a href="theory/index.html">Theory manual</a> when you want the mathematical meaning behind the vectorized implementation.</p>
    </div>
  </div>
  <div class="fp-step">
    <span class="fp-step-num">4</span>
    <div class="fp-step-content">
      <h3>Look Up Exact Behavior</h3>
      <p>Use the <a href="api.html">API reference</a> when you need precise call signatures, return types, and module ownership.</p>
    </div>
  </div>
</div>

<div class="fp-callout-grid">
  <div class="fp-callout fp-callout--blue">
    <div class="fp-callout-icon">
      <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M22 12h-4l-3 9L9 3l-3 9H2"/></svg>
    </div>
    <div class="fp-callout-body">
      <strong>Three-layer architecture</strong>
      <p>The <a href="manual/index.html">Guide</a> teaches workflows. The <a href="theory/index.html">Theory</a> explains the math. The <a href="api.html">API</a> gives exact signatures.</p>
    </div>
  </div>
  <div class="fp-callout fp-callout--amber">
    <div class="fp-callout-icon">
      <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M10.29 3.86L1.82 18a2 2 0 001.71 3h16.94a2 2 0 001.71-3L13.71 3.86a2 2 0 00-3.42 0z"/><line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/></svg>
    </div>
    <div class="fp-callout-body">
      <strong>One-based vs zero-based</strong>
      <p>Public tables use <strong>1-based</strong> node numbering (FemLab convention). Internal NumPy arrays are 0-based. Both manuals explain where the shift happens.</p>
    </div>
  </div>
</div>
