API Reference
=============

.. raw:: html

   <div class="fp-hero fp-hero--compact">
     <div class="fp-hero-title">API Reference</div>
     <div class="fp-hero-subtitle">
       Module-level documentation &mdash; find symbols, understand ownership, read full signatures.
     </div>
   </div>

.. container:: fp-card-grid

   .. container:: fp-card fp-card--accent-blue

      **Core Workflow**

      - ``femlabpy.core`` &mdash; global array allocation
      - ``femlabpy.assembly`` &mdash; element scatter
      - ``femlabpy.loads`` &mdash; nodal force mapping
      - ``femlabpy.boundary`` &mdash; prescribed displacements

   .. container:: fp-card fp-card--accent-red

      **Dynamics & Nonlinear**

      - ``femlabpy.dynamics`` &mdash; time-history solvers, FRF
      - ``femlabpy.damping`` &mdash; Rayleigh & modal damping
      - ``femlabpy.modal`` &mdash; eigenproblem, participation
      - ``femlabpy.solvers`` &mdash; legacy nonlinear drivers

   .. container:: fp-card fp-card--accent-green

      **Element Kernels**

      - ``femlabpy.elements.bars`` &mdash; bar K, f, M
      - ``femlabpy.elements.triangles`` &mdash; T3 family
      - ``femlabpy.elements.quads`` &mdash; Q4 elastic/plastic/mass
      - ``femlabpy.elements.solids`` &mdash; T4 & H8 solids

   .. container:: fp-card

      **Support & I/O**

      - ``femlabpy.periodic`` &mdash; pairing & homogenization
      - ``femlabpy.postprocess`` &mdash; reactions, result tables
      - ``femlabpy.plotting`` &mdash; Matplotlib views
      - ``femlabpy.io.gmsh`` &mdash; Gmsh mesh import

.. raw:: html

   <div class="fp-callout-grid">
     <div class="fp-callout fp-callout--blue">
       <div class="fp-callout-icon">
         <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M2 3h6a4 4 0 014 4v14a3 3 0 00-3-3H2z"/><path d="M22 3h-6a4 4 0 00-4 4v14a3 3 0 013-3h7z"/></svg>
       </div>
       <div class="fp-callout-body">
         <strong>Reading order</strong>
         <p>Read by responsibility: <code>core</code> &rarr; <code>assembly / loads / boundary</code> &rarr; element kernels &rarr; dynamics &rarr; periodic / I/O.</p>
       </div>
     </div>
   </div>

Complete Module Index
---------------------

.. autosummary::
   :toctree: generated/

   femlabpy.core
   femlabpy.assembly
   femlabpy.boundary
   femlabpy.loads
   femlabpy.solvers
   femlabpy.dynamics
   femlabpy.damping
   femlabpy.modal
   femlabpy.periodic
   femlabpy.materials
   femlabpy.elements.bars
   femlabpy.elements.quads
   femlabpy.elements.triangles
   femlabpy.elements.solids
   femlabpy.postprocess
   femlabpy.plotting
   femlabpy.matlab
   femlabpy.io.gmsh
