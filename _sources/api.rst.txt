API Reference
=============

.. container:: fp-home-intro

   .. container:: fp-kicker

      Module-level ownership and full symbol reference

   .. container:: fp-lead

      This section is the module-level reference for ``femlabpy``. It is
      organized to help with two tasks at the same time: finding the exact
      symbol you need and understanding which module owns which part of the
      workflow.

.. container:: fp-card-grid

   .. container:: fp-card

      **Core workflow modules**

      - ``femlabpy.core`` allocates the global containers.
      - ``femlabpy.assembly`` scatters element matrices and vectors.
      - ``femlabpy.loads`` maps nodal force tables into the global right-hand side.
      - ``femlabpy.boundary`` applies prescribed displacements and exact constraints.

   .. container:: fp-card

      **Dynamics and nonlinear modules**

      - ``femlabpy.dynamics`` contains load builders, time-history solvers, FRF
        utilities, and response plotting helpers.
      - ``femlabpy.damping`` builds Rayleigh and modal damping matrices.
      - ``femlabpy.modal`` solves the eigenproblem and computes participation.
      - ``femlabpy.solvers`` exposes the larger legacy nonlinear drivers.

   .. container:: fp-card

      **Element kernel modules**

      - ``femlabpy.elements.bars`` implements bar stiffness, internal-force, and mass routines.
      - ``femlabpy.elements.triangles`` implements the T3 family.
      - ``femlabpy.elements.quads`` implements Q4 elastic, scalar, plastic, and mass kernels.
      - ``femlabpy.elements.solids`` implements T4 and H8 solid kernels.

   .. container:: fp-card

      **Support modules**

      - ``femlabpy.periodic`` handles periodic pairing, constraints, and homogenization.
      - ``femlabpy.postprocess`` extracts reactions and compact result tables.
      - ``femlabpy.plotting`` provides lightweight Matplotlib views.
      - ``femlabpy.io.gmsh`` reads Gmsh meshes into normalized tables.
      - ``femlabpy.matlab`` exposes packaged legacy benchmark decks and wrappers.
      - ``femlabpy.materials`` re-exports invariant and plasticity helpers.

How To Read This Section
------------------------

If you are new to the codebase, do not read the module pages in alphabetical
order. Read them by responsibility:

1. ``femlabpy.core`` for global array allocation.
2. ``femlabpy.assembly``, ``femlabpy.loads``, and ``femlabpy.boundary`` for
   system construction.
3. The element modules for low-level stiffness, internal-force, and mass kernels.
4. ``femlabpy.modal``, ``femlabpy.damping``, ``femlabpy.dynamics``, and
   ``femlabpy.solvers`` for modal, transient, and nonlinear workflows.
5. ``femlabpy.periodic`` and ``femlabpy.io.gmsh`` for periodic or mesh-driven
   workflows.

.. container:: fp-section-note

   Each module page now includes both a summary table and the full rendered
   documentation for its public functions and classes. The summary table is the
   fastest way to find a symbol, and the lower part of each page holds the full
   reference text.

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
