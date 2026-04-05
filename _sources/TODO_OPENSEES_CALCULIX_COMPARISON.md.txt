# femlabpy vs OpenSees vs CalculiX: Comprehensive Gap Analysis and Implementation Roadmap

**Document Version:** 1.0  
**Created:** 2026-04-05  
**Status:** Active Development Roadmap  
**Target Completion:** Long-term (12-24 months)  

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Project Overview](#2-project-overview)
3. [OpenSees Source Code Analysis](#3-opensees-source-code-analysis)
   - 3.1 [Repository Structure](#31-repository-structure)
   - 3.2 [Element Library Analysis](#32-element-library-analysis)
   - 3.3 [Material Library Analysis](#33-material-library-analysis)
   - 3.4 [Analysis Framework Analysis](#34-analysis-framework-analysis)
   - 3.5 [Integration Algorithms](#35-integration-algorithms)
4. [CalculiX Source Code Analysis](#4-calculix-source-code-analysis)
   - 4.1 [Solver Architecture](#41-solver-architecture)
   - 4.2 [Element Types](#42-element-types)
   - 4.3 [Material Models](#43-material-models)
   - 4.4 [Analysis Capabilities](#44-analysis-capabilities)
5. [femlabpy Current State Analysis](#5-femlabpy-current-state-analysis)
   - 5.1 [Existing Elements](#51-existing-elements)
   - 5.2 [Existing Materials](#52-existing-materials)
   - 5.3 [Existing Solvers](#53-existing-solvers)
   - 5.4 [Existing I/O](#54-existing-io)
6. [Gap Analysis Matrix](#6-gap-analysis-matrix)
7. [Implementation Roadmap](#7-implementation-roadmap)
   - Phase 1: Core Structural Elements
   - Phase 2: Material Library Expansion
   - Phase 3: Advanced Analysis
   - Phase 4: High-Performance Computing
   - Phase 5: Advanced Elements
8. [Detailed Implementation Tasks](#8-detailed-implementation-tasks)
9. [Testing and Validation Strategy](#9-testing-and-validation-strategy)
10. [References](#10-references)

---

## 1. Executive Summary

This document provides an exhaustive analysis comparing `femlabpy` (a Python finite element library for teaching) against two major open-source FEM packages:

1. **OpenSees** (Open System for Earthquake Engineering Simulation) - Berkeley, CA
2. **CalculiX** - A free 3D structural FEM program by Guido Dhondt and Klaus Wittig

The analysis identifies **247 specific implementation tasks** organized into **5 development phases** spanning approximately **18-24 months** of development effort.

### Key Findings

| Category | OpenSees | CalculiX | femlabpy | Gap |
|----------|----------|----------|----------|-----|
| Element Types | 47 categories | 40+ types | 6 types | ~85% missing |
| Material Models | 170+ uniaxial | 30+ models | 2 models | ~95% missing |
| Analysis Types | 12+ types | 15+ types | 5 types | ~70% missing |
| Integration Schemes | 15+ | 10+ | 4 | ~70% missing |

---

## 2. Project Overview

### 2.1 femlabpy Mission Statement

`femlabpy` is a Python FEM library designed for:
- Teaching finite element concepts
- Providing clean, readable implementations
- Matching classroom MATLAB/Scilab FemLab workflows
- Enabling research prototyping

### 2.2 OpenSees Mission Statement

OpenSees is designed for:
- Earthquake engineering simulation
- Nonlinear structural analysis
- Soil-structure interaction
- Performance-based earthquake engineering
- Research and professional practice

### 2.3 CalculiX Mission Statement

CalculiX is designed for:
- 3D structural analysis
- Thermal analysis
- CFD pre/post-processing compatibility
- Industrial applications (aerospace, automotive)
- Linear and nonlinear analysis

### 2.4 Development Philosophy Comparison

| Aspect | femlabpy | OpenSees | CalculiX |
|--------|----------|----------|----------|
| Language | Python | C++/Tcl/Python | Fortran/C |
| Target User | Students/Researchers | Researchers/Engineers | Engineers/Industry |
| Code Style | Readable/Educational | Production | Production |
| Performance | NumPy/SciPy | Optimized C++ | Fortran optimized |
| Documentation | In-code | Wiki/Manual | PDF Manual |

---

## 3. OpenSees Source Code Analysis

### 3.1 Repository Structure

**Repository:** https://github.com/OpenSees/OpenSees

```
OpenSees/
├── SRC/
│   ├── actor/              # Parallel computing actors
│   ├── analysis/           # Analysis framework
│   │   ├── algorithm/      # Solution algorithms (Newton, BFGS, etc.)
│   │   ├── analysis/       # Analysis types
│   │   ├── dof_grp/        # DOF grouping
│   │   ├── fe_ele/         # FE element wrappers
│   │   ├── handler/        # Constraint handlers
│   │   ├── integrator/     # Time integration schemes
│   │   ├── model/          # Analysis model
│   │   └── numberer/       # DOF numbering
│   ├── convergenceTest/    # Convergence criteria
│   ├── coordTransformation/# Coordinate transformations
│   ├── damage/             # Damage models
│   ├── damping/            # Damping models
│   ├── domain/             # Domain management
│   ├── element/            # Element library (47 categories!)
│   ├── graph/              # Graph algorithms
│   ├── material/           # Material library
│   │   ├── nD/             # Multi-dimensional materials
│   │   ├── section/        # Section definitions
│   │   ├── uniaxial/       # 1D materials (170+ types!)
│   │   └── yieldSurface/   # Yield surfaces
│   ├── matrix/             # Matrix operations
│   ├── recorder/           # Output recorders
│   ├── reliability/        # Reliability analysis
│   ├── system_of_eqn/      # Equation solvers
│   └── utility/            # Utilities
├── EXAMPLES/               # Example models
└── tests/                  # Test suite
```

### 3.2 Element Library Analysis

OpenSees contains **47 element categories** with hundreds of specific implementations:

#### 3.2.1 Truss and Bar Elements
```
SRC/element/truss/
├── Truss.cpp                    # Basic truss element
├── TrussSection.cpp             # Truss with section
├── CorotTruss.cpp               # Corotational truss
└── CorotTrussSection.cpp        # Corotational with section
```

**femlabpy equivalent:** `elements/bars.py` - Basic truss only
**Gap:** Corotational formulation, section-based truss

#### 3.2.2 Beam Elements (2D)
```
SRC/element/beam2d/
├── beam2d02.cpp                 # 2D Euler-Bernoulli beam
├── beam2d03.cpp                 # 2D beam with P-Delta
└── beam2d04.cpp                 # 2D beam with distributed loads
```

**femlabpy equivalent:** NONE
**Gap:** Complete 2D beam element family missing

#### 3.2.3 Beam Elements (3D)
```
SRC/element/beam3d/
├── beam3d01.cpp                 # Basic 3D beam
└── beam3d02.cpp                 # 3D beam with P-Delta
```

**femlabpy equivalent:** NONE
**Gap:** Complete 3D beam element family missing

#### 3.2.4 Elastic Beam-Column Elements
```
SRC/element/elasticBeamColumn/
├── ElasticBeam2d.cpp            # 2D elastic beam-column
├── ElasticBeam3d.cpp            # 3D elastic beam-column
├── ElasticTimoshenkoBeam2d.cpp  # 2D Timoshenko beam
├── ElasticTimoshenkoBeam3d.cpp  # 3D Timoshenko beam
├── ModElasticBeam2d.cpp         # Modified elastic beam
└── ModElasticBeam3d.cpp         # Modified 3D beam
```

**femlabpy equivalent:** NONE
**Gap:** Entire elastic beam-column family missing

#### 3.2.5 Displacement-Based Beam-Column
```
SRC/element/dispBeamColumn/
├── DispBeamColumn2d.cpp         # 2D displacement-based beam
├── DispBeamColumn3d.cpp         # 3D displacement-based beam
├── DispBeamColumn2dThermal.cpp  # Thermal coupling
├── DispBeamColumn3dThermal.cpp  # 3D thermal
├── DispBeamColumnNL2d.cpp       # Nonlinear 2D
└── TimoshenkoBeamColumn2d.cpp   # Timoshenko formulation
```

**femlabpy equivalent:** NONE
**Gap:** Displacement-based beam formulation missing

#### 3.2.6 Force-Based Beam-Column
```
SRC/element/forceBeamColumn/
├── ForceBeamColumn2d.cpp        # 2D force-based beam
├── ForceBeamColumn3d.cpp        # 3D force-based beam
├── ForceBeamColumnWarping2d.cpp # With warping
├── ForceBeamColumnCBDI2d.cpp    # CBDI formulation
└── ElasticForceBeamColumn2d.cpp # Elastic version
```

**femlabpy equivalent:** NONE
**Gap:** Force-based beam formulation missing

#### 3.2.7 Beam with Hinges
```
SRC/element/beamWithHinges/
├── BeamWithHinges2d.cpp         # 2D beam with plastic hinges
├── BeamWithHinges3d.cpp         # 3D beam with plastic hinges
└── TclBeamWithHingesCommand.cpp # Tcl interface
```

**femlabpy equivalent:** NONE
**Gap:** Plastic hinge beam elements missing

#### 3.2.8 Quadrilateral Elements
```
SRC/element/fourNodeQuad/
├── ConstantPressureVolumeQuad.cpp    # Incompressible
├── EnhancedQuad.cpp                   # Enhanced strain
├── FourNodeQuad.cpp                   # Standard Q4
├── FourNodeQuad3d.cpp                 # 3D membrane
├── FourNodeQuadUP.cpp                 # Undrained (u-p)
├── FourNodeQuadWithSensitivity.cpp    # With sensitivity
├── NineNodeMixedQuad.cpp              # 9-node mixed
├── NineNodeQuad.cpp                   # 9-node (Q9)
├── SixNodeTri.cpp                     # 6-node triangle (T6)
└── SSPquad.cpp                        # Stabilized single-point
```

**femlabpy equivalent:** `elements/quads.py` - Basic Q4 only
**Gap:** Q9, T6, enhanced strain, u-p formulation, mixed methods

#### 3.2.9 Triangle Elements
```
SRC/element/triangle/
├── Tri31.cpp                    # 3-node triangle
├── SixNodeTri.cpp               # 6-node triangle
└── TripleFrictionPendulum.cpp   # Special friction element
```

**femlabpy equivalent:** `elements/triangles.py` - T3 only
**Gap:** T6, special triangles

#### 3.2.10 Brick (Hexahedron) Elements
```
SRC/element/brick/
├── Brick.cpp                    # 8-node brick (H8)
├── BbarBrick.cpp                # B-bar formulation
├── BrickUP.cpp                  # Undrained u-p
├── Brick20N.cpp                 # 20-node brick (H20)
├── TwentySevenNodeBrick.cpp     # 27-node brick (H27)
├── TwentyEightNodeBrickUP.cpp   # 28-node u-p
├── SSPbrick.cpp                 # Stabilized single-point
└── SSPbrickUP.cpp               # SSP undrained
```

**femlabpy equivalent:** `elements/solids.py` - Basic H8 only
**Gap:** H20, H27, B-bar, u-p formulation, SSP

#### 3.2.11 Tetrahedron Elements
```
SRC/element/tetrahedron/
├── FourNodeTetrahedron.cpp      # 4-node tet (T4)
├── TenNodeTetrahedron.cpp       # 10-node tet (T10)
└── TenNodeTetrahedronThermal.cpp # Thermal T10
```

**femlabpy equivalent:** `elements/solids.py` - Basic T4 only
**Gap:** T10, thermal coupling

#### 3.2.12 Shell Elements
```
SRC/element/shell/
├── ASDShellQ4.cpp               # MITC4 shell
├── ShellANDeS.cpp               # ANDES formulation
├── ShellDKGQ.cpp                # DKG quadrilateral
├── ShellDKGT.cpp                # DKG triangle
├── ShellMITC4.cpp               # MITC4 shell
├── ShellMITC4Thermal.cpp        # Thermal MITC4
├── ShellMITC9.cpp               # MITC9 shell
├── ShellNL.cpp                  # Nonlinear shell
├── ShellNLDKGQ.cpp              # Nonlinear DKGQ
└── ShellNLDKGT.cpp              # Nonlinear DKGT
```

**femlabpy equivalent:** NONE
**Gap:** Complete shell element family missing

#### 3.2.13 Zero-Length Elements
```
SRC/element/zeroLength/
├── ZeroLength.cpp               # Zero-length spring
├── ZeroLengthContact2D.cpp      # 2D contact
├── ZeroLengthContact3D.cpp      # 3D contact
├── ZeroLengthContactASDimplex.cpp # ASD contact
├── ZeroLengthContactNTS2D.cpp   # NTS contact
├── ZeroLengthND.cpp             # Multi-dimensional
├── ZeroLengthRocking.cpp        # Rocking element
└── ZeroLengthSection.cpp        # Section-based
```

**femlabpy equivalent:** NONE
**Gap:** Complete zero-length family missing (springs, contacts)

#### 3.2.14 Joint Elements
```
SRC/element/joint/
├── BeamColumnJoint2d.cpp        # 2D beam-column joint
├── BeamColumnJoint3d.cpp        # 3D beam-column joint
├── ElasticTubularJoint.cpp      # Tubular joint
├── Joint2D.cpp                  # 2D joint
├── Joint3D.cpp                  # 3D joint
└── MP_Joint2D.cpp               # Multi-point joint
```

**femlabpy equivalent:** NONE
**Gap:** Complete joint element family missing

#### 3.2.15 Bearing Elements
```
SRC/element/elastomericBearing/
├── ElastomericBearingBoucWen2d.cpp   # Bouc-Wen 2D
├── ElastomericBearingBoucWen3d.cpp   # Bouc-Wen 3D
├── ElastomericBearingPlasticity2d.cpp # Plasticity 2D
├── ElastomericBearingPlasticity3d.cpp # Plasticity 3D
├── ElastomericBearingUFRP2d.cpp      # UFRP 2D
├── ElastomericX.cpp                   # ElastomericX
└── HDR.cpp                            # High damping rubber

SRC/element/frictionBearing/
├── FlatSliderSimple2d.cpp       # Flat slider 2D
├── FlatSliderSimple3d.cpp       # Flat slider 3D
├── FPBearingPTV.cpp             # FP bearing
├── SingleFPSimple2d.cpp         # Single FP 2D
├── SingleFPSimple3d.cpp         # Single FP 3D
├── RJWatsonEqsBearing.cpp       # RJ Watson bearing
├── TripleFrictionPendulum.cpp   # Triple FP
└── TripleFrictionPendulumX.cpp  # Triple FPX
```

**femlabpy equivalent:** NONE
**Gap:** Complete bearing element families missing

#### 3.2.16 Cable Elements
```
SRC/element/catenaryCable/
├── CatenaryCable.cpp            # Catenary cable
└── CorotationalCable.cpp        # Corotational cable
```

**femlabpy equivalent:** NONE
**Gap:** Cable elements missing

#### 3.2.17 Link Elements
```
SRC/element/twoNodeLink/
├── TwoNodeLink.cpp              # General two-node link
├── TwoNodeLinkSection.cpp       # Section-based link
└── Inerter.cpp                  # Inerter element
```

**femlabpy equivalent:** NONE
**Gap:** Link elements missing

#### 3.2.18 MVLEM (Multiple-Vertical-Line-Element-Model)
```
SRC/element/mvlem/
├── MVLEM.cpp                    # Basic MVLEM
├── MVLEM_3D.cpp                 # 3D MVLEM
├── SFI_MVLEM.cpp                # Shear-flexure interaction
└── SFI_MVLEM_3D.cpp             # 3D SFI
```

**femlabpy equivalent:** NONE
**Gap:** Wall modeling elements missing

#### 3.2.19 Masonry Elements
```
SRC/element/masonry/
├── BeamGT.cpp                   # Masonry beam
└── MasonPan12.cpp               # Masonry panel
```

**femlabpy equivalent:** NONE
**Gap:** Masonry elements missing

#### 3.2.20 Pipe Elements
```
SRC/element/pipe/
├── Pipe.cpp                     # Pipe element
└── ElasticPipe.cpp              # Elastic pipe
```

**femlabpy equivalent:** NONE
**Gap:** Pipe elements missing

#### 3.2.21 Absorbent Boundaries
```
SRC/element/absorbentBoundaries/
├── ASDAbsorbingBoundary2D.cpp   # 2D absorbing BC
└── ASDAbsorbingBoundary3D.cpp   # 3D absorbing BC
```

**femlabpy equivalent:** NONE
**Gap:** Absorbing boundaries missing

#### 3.2.22 PML (Perfectly Matched Layer)
```
SRC/element/PML/
├── PML3D.cpp                    # 3D PML
├── PML2D.cpp                    # 2D PML
└── PML2DVISCOUS.cpp             # Viscous PML
```

**femlabpy equivalent:** NONE
**Gap:** PML elements missing

#### 3.2.23 Surface Load Elements
```
SRC/element/surfaceLoad/
├── BrickSurfaceLoad.cpp         # Brick surface
├── QuadSurfaceLoad.cpp          # Quad surface
├── SurfaceLoad.cpp              # Generic surface
├── TriSurfaceLoad.cpp           # Triangle surface
└── VS3D4.cpp                    # Viscous surface
```

**femlabpy equivalent:** NONE
**Gap:** Surface load elements missing

#### 3.2.24 Updated Lagrangian Beam-Column
```
SRC/element/updatedLagrangianBeamColumn/
├── Beam2dPointLoad.cpp          # Point load beam
├── Beam2dTempLoad.cpp           # Temperature load
├── Beam2dUniformLoad.cpp        # Uniform load
├── Elastic2DGNL.cpp             # Geometric NL
├── InelasticYS2DGNL.cpp         # Inelastic yield surface
├── UpdatedLagrangianBeam2D.cpp  # 2D updated Lagrangian
└── UpdatedLagrangianBeam3D.cpp  # 3D updated Lagrangian
```

**femlabpy equivalent:** NONE
**Gap:** Updated Lagrangian formulation missing

#### 3.2.25 Gradient Inelastic Beam-Column
```
SRC/element/gradientInelasticBeamColumn/
├── GradientInelasticBeamColumn2d.cpp
└── GradientInelasticBeamColumn3d.cpp
```

**femlabpy equivalent:** NONE
**Gap:** Gradient-enhanced plasticity missing

#### 3.2.26 Mixed Beam-Column
```
SRC/element/mixedBeamColumn/
├── MixedBeamColumn2d.cpp        # Mixed 2D
├── MixedBeamColumn3d.cpp        # Mixed 3D
└── MixedBeamColumnAsym3d.cpp    # Asymmetric 3D
```

**femlabpy equivalent:** NONE
**Gap:** Mixed formulation missing

### 3.3 Material Library Analysis

OpenSees contains **170+ uniaxial materials** and extensive nD materials:

#### 3.3.1 Steel Materials (Uniaxial)
```
SRC/material/uniaxial/
├── Steel01.cpp                  # Bilinear steel with kinematic hardening
├── Steel01Thermal.cpp           # Thermal version
├── Steel02.cpp                  # Giuffré-Menegotto-Pinto steel
├── Steel02Fatigue.cpp           # With fatigue
├── Steel02Thermal.cpp           # Thermal version
├── Steel03.cpp                  # With asymmetric properties
├── Steel4.cpp                   # Advanced steel model
├── SteelBRB.cpp                 # Buckling-restrained brace
├── SteelDRC.cpp                 # Damage-rate-dependent
├── SteelECThermal.cpp           # Eurocode thermal
├── SteelFractureDI.cpp          # Fracture with damage index
├── SteelMP.cpp                  # Menegotto-Pinto
├── SteelMPF.cpp                 # MPF steel
├── ReinforcingSteel.cpp         # Reinforcing bars
├── DoddRestr.cpp                # Dodd-Restrepo model
├── Dodd_Restrepo.cpp            # Alternative implementation
└── PrestressedSteelMaterial.cpp # Prestressing steel
```

**femlabpy equivalent:** NONE
**Gap:** Complete steel material library missing

#### 3.3.2 Concrete Materials (Uniaxial)
```
SRC/material/uniaxial/
├── Concrete01.cpp               # Kent-Scott-Park unconfined
├── Concrete01WithSITC.cpp       # With strain-induced tension
├── Concrete02.cpp               # Linear tension softening
├── Concrete02IS.cpp             # Initial stress
├── Concrete02Thermal.cpp        # Thermal version
├── Concrete04.cpp               # Popovics confined/unconfined
├── Concrete05.cpp               # Confined concrete
├── Concrete06.cpp               # Advanced model
├── Concrete07.cpp               # Chang & Mander
├── ConcreteCM.cpp               # Confined Model
├── ConcreteD.cpp                # Damage-based
├── ConcreteECThermal.cpp        # Eurocode thermal
├── ConcreteSakaiKawashima.cpp   # Sakai-Kawashima model
├── ConcreteZBH_fitted.cpp       # ZBH fitted
├── ConcreteZBH_original.cpp     # ZBH original
├── ConcreteZBH_smoothed.cpp     # ZBH smoothed
├── ConcretewBeta.cpp            # With beta parameter
├── ConfinedConcrete01.cpp       # Legeron-Paultre model
├── ECC01.cpp                    # Engineered cementitious
├── FRCC.cpp                     # Fiber-reinforced concrete
├── FRPConfinedConcrete.cpp      # FRP confined
├── FRPConfinedConcrete02.cpp    # Alternative FRP
├── TDConcrete.cpp               # Time-dependent
├── TDConcreteEXP.cpp            # Experimental TD
├── TDConcreteMC10.cpp           # Model Code 2010
├── TDConcreteMC10NL.cpp         # MC10 nonlinear
└── TDConcreteNL.cpp             # Nonlinear TD
```

**femlabpy equivalent:** NONE
**Gap:** Complete concrete material library missing

#### 3.3.3 Hysteretic Materials
```
SRC/material/uniaxial/
├── BoucWenMaterial.cpp          # Bouc-Wen model
├── BoucWenOriginal.cpp          # Original Bouc-Wen
├── BoucWenInfill.cpp            # Infill wall
├── BWBN.cpp                     # Bouc-Wen-Baber-Noori
├── HystereticMaterial.cpp       # Generic hysteretic
├── HystereticAsym.cpp           # Asymmetric hysteretic
├── HystereticPoly.cpp           # Polynomial hysteretic
├── HystereticSmooth.cpp         # Smooth hysteretic
├── HystereticSMMaterial.cpp     # Smooth multilinear
├── OOHystereticMaterial.cpp     # Object-oriented hysteretic
├── Pinching4Material.cpp        # Pinching with 4 parameters
├── DegradingPinchedBW.cpp       # Degrading pinched
├── ModIMKPeakOriented.cpp       # Modified IMK peak-oriented
├── ModIMKPeakOriented02.cpp     # IMK version 2
├── ModIMKPinching.cpp           # IMK pinching
├── ModIMKPinching02.cpp         # IMK pinching v2
├── IMKBilin.cpp                 # IMK bilinear
├── IMKPeakOriented.cpp          # IMK peak-oriented
├── IMKPinching.cpp              # IMK pinching
├── Bilin.cpp                    # Bilinear with stiffness degradation
├── Bilin02.cpp                  # Bilin version 2
├── ShearPanelMaterial.cpp       # Shear panel
└── SAWSMaterial.cpp             # SAWS (wood)
```

**femlabpy equivalent:** NONE
**Gap:** Complete hysteretic material library missing

#### 3.3.4 Damper Materials
```
SRC/material/uniaxial/
├── ViscousDamper.cpp            # Viscous damper
├── ViscousMaterial.cpp          # Viscous material
├── BilinearOilDamper.cpp        # Oil damper
├── DamperMaterial.cpp           # Generic damper
├── Maxwell.cpp                  # Maxwell model
├── ViscoelasticGap.cpp          # Viscoelastic gap
├── CoulombDamperMaterial.cpp    # Coulomb damper
├── Hertzdamp.cpp                # Hertz damping
├── APDFMD.cpp                   # Adaptive passive
├── APDMD.cpp                    # APD material damper
└── APDVFD.cpp                   # APD viscous fluid
```

**femlabpy equivalent:** NONE
**Gap:** Complete damper material library missing

#### 3.3.5 Gap and Contact Materials
```
SRC/material/uniaxial/
├── ElasticPPMaterial.cpp        # Elastic-perfectly-plastic
├── EPPGapMaterial.cpp           # EPP with gap
├── HookGap.cpp                  # Hook with gap
├── HyperbolicGapMaterial.cpp    # Hyperbolic gap
├── ImpactMaterial.cpp           # Impact material
├── JankowskiImpact.cpp          # Jankowski impact
└── Ratchet.cpp                  # Ratchet material
```

**femlabpy equivalent:** NONE
**Gap:** Gap and contact materials missing

#### 3.3.6 SMA (Shape Memory Alloy) Materials
```
SRC/material/uniaxial/
├── SMAMaterial.cpp              # Shape memory alloy
├── SelfCenteringMaterial.cpp    # Self-centering
├── ASD_SMA_3K.cpp               # ASD SMA 3K
└── FlagShapeMaterial.cpp        # Flag-shaped hysteresis
```

**femlabpy equivalent:** NONE
**Gap:** SMA materials missing

#### 3.3.7 Isolator Materials
```
SRC/material/uniaxial/
├── KikuchiAikenHDR.cpp          # HDR isolator
├── KikuchiAikenLRB.cpp          # Lead rubber bearing
└── Neoprene.cpp                 # Neoprene material
```

**femlabpy equivalent:** NONE
**Gap:** Isolator materials missing

#### 3.3.8 Wrapper/Modifier Materials
```
SRC/material/uniaxial/
├── ParallelMaterial.cpp         # Parallel combination
├── SeriesMaterial.cpp           # Series combination
├── MultiplierMaterial.cpp       # Multiplier wrapper
├── MinMaxMaterial.cpp           # Min-max limits
├── InitStrainMaterial.cpp       # Initial strain
├── InitStressMaterial.cpp       # Initial stress
├── FatigueMaterial.cpp          # Fatigue wrapper
├── TensionOnlyMaterial.cpp      # Tension only
├── SimpleFractureMaterial.cpp   # Simple fracture
├── DuctileFracture.cpp          # Ductile fracture
├── PathIndependentMaterial.cpp  # Path independent
├── PenaltyMaterial.cpp          # Penalty material
└── BackboneMaterial.cpp         # Backbone curve
```

**femlabpy equivalent:** NONE
**Gap:** Material wrappers missing

#### 3.3.9 Soil p-y, t-z, Q-z Materials
```
SRC/material/uniaxial/
├── QbSandCPT.cpp                # Q-b sand CPT
├── TzSandCPT.cpp                # T-z sand CPT
└── pyUCLA.cpp                   # UCLA p-y curves
```

**femlabpy equivalent:** NONE
**Gap:** Soil-structure interaction materials missing

#### 3.3.10 Multi-Dimensional (nD) Materials
```
SRC/material/nD/
├── ElasticIsotropicMaterial.cpp        # Elastic isotropic
├── ElasticIsotropic3D.cpp              # 3D elastic
├── ElasticIsotropicPlaneStrain2D.cpp   # Plane strain
├── ElasticIsotropicPlaneStress2D.cpp   # Plane stress
├── ElasticOrthotropicMaterial.cpp      # Orthotropic
├── J2Plasticity.cpp                    # J2 plasticity
├── PlasticDamageConcrete3d.cpp         # Plastic damage concrete
├── PlasticDamageConcretePlaneStress.cpp # Plane stress PDC
├── DruckerPrager.cpp                   # Drucker-Prager
├── DruckerPrager3D.cpp                 # 3D DP
├── DruckerPragerPlaneStrain.cpp        # Plane strain DP
├── CycLiqCP.cpp                        # Cyclic liquefaction
├── CycLiqCP3D.cpp                      # 3D cyclic liquefaction
├── CycLiqCPSP.cpp                      # SP version
├── CycLiqCPSP3D.cpp                    # 3D SP
├── PressureDependMultiYield.cpp        # PDMY soil
├── PressureDependMultiYield02.cpp      # PDMY02
├── PressureIndependMultiYield.cpp      # PIMY soil
├── ContactMaterial2D.cpp               # 2D contact
├── ContactMaterial3D.cpp               # 3D contact
├── InitialStateAnalysisWrapper.cpp     # Initial state
├── ManzariDafalias.cpp                 # Manzari-Dafalias sand
├── ManzariDafalias3D.cpp               # 3D MD
├── ManzariDafaliasPlaneStrain.cpp      # PS MD
├── PM4Sand.cpp                         # PM4Sand model
├── PM4Silt.cpp                         # PM4Silt model
├── StressDensityModel.cpp              # Stress density
└── J2BeamFiber2d.cpp                   # Beam fiber J2
```

**femlabpy equivalent:** `materials/plasticity.py` - von Mises only
**Gap:** ~95% of nD material library missing

#### 3.3.11 Section Models
```
SRC/material/section/
├── ElasticSection2d.cpp         # 2D elastic section
├── ElasticSection3d.cpp         # 3D elastic section
├── ElasticShearSection2d.cpp    # With shear 2D
├── ElasticShearSection3d.cpp    # With shear 3D
├── ElasticTubeSection3d.cpp     # Tube section
├── ElasticWarpingShearSection2d.cpp  # Warping section
├── FiberSection2d.cpp           # 2D fiber section
├── FiberSection3d.cpp           # 3D fiber section
├── FiberSectionAsym3d.cpp       # Asymmetric fiber
├── FiberSectionGJ.cpp           # GJ fiber
├── FiberSectionGJThermal.cpp    # Thermal fiber
├── FiberSectionShear3d.cpp      # Shear fiber
├── FiberSectionWarping3d.cpp    # Warping fiber
├── GenericSection1d.cpp         # Generic 1D
├── GenericSectionNd.cpp         # Generic nD
├── LayeredShellFiberSection.cpp # Layered shell
├── LayeredShellFiberSectionThermal.cpp # Thermal layered
├── MembranePlateFiberSection.cpp # Membrane plate
├── NDFiberSection2d.cpp         # nD fiber 2D
├── NDFiberSection3d.cpp         # nD fiber 3D
├── NDFiberSectionWarping2d.cpp  # nD warping 2D
├── ParallelSection.cpp          # Parallel combination
├── SectionAggregator.cpp        # Section aggregation
├── UniaxialSection2d.cpp        # Uniaxial 2D
└── WideFlangeSectionIntegration.cpp # Wide flange
```

**femlabpy equivalent:** NONE
**Gap:** Complete section model library missing

### 3.4 Analysis Framework Analysis

#### 3.4.1 Analysis Types
```
SRC/analysis/analysis/
├── DirectIntegrationAnalysis.cpp     # Time history
├── EigenAnalysis.cpp                 # Eigenvalue
├── DomainDecompositionAnalysis.cpp   # Parallel
├── StaticAnalysis.cpp                # Static
├── StaticDomainDecompositionAnalysis.cpp # Parallel static
├── SubstructuringAnalysis.cpp        # Substructuring
├── TransientAnalysis.cpp             # Transient
├── TransientDomainDecompositionAnalysis.cpp # Parallel transient
├── VariableTimeStepDirectIntegrationAnalysis.cpp # Variable dt
└── PFEMAnalysis.cpp                  # PFEM
```

**femlabpy equivalent:** Basic static and modal
**Gap:** Domain decomposition, PFEM, variable timestep

#### 3.4.2 Solution Algorithms
```
SRC/analysis/algorithm/
├── Linear.cpp                   # Linear algorithm
├── NewtonRaphson.cpp            # Newton-Raphson
├── NewtonLineSearch.cpp         # Newton with line search
├── ModifiedNewton.cpp           # Modified Newton
├── BFGS.cpp                     # BFGS quasi-Newton
├── Broyden.cpp                  # Broyden quasi-Newton
├── KrylovNewton.cpp             # Krylov-Newton
├── PeriodicNewton.cpp           # Periodic Newton
├── SecantNewton.cpp             # Secant Newton
├── RegulaFalsiLineSearch.cpp    # Regula Falsi
├── BisectionLineSearch.cpp      # Bisection line search
└── InitialInterpolatedLineSearch.cpp # Interpolated LS
```

**femlabpy equivalent:** Basic Newton-Raphson
**Gap:** Line search, quasi-Newton methods, Krylov

#### 3.4.3 Constraint Handlers
```
SRC/analysis/handler/
├── ConstraintHandler.cpp        # Base class
├── LagrangeConstraintHandler.cpp # Lagrange multipliers
├── PenaltyConstraintHandler.cpp # Penalty method
├── PlainHandler.cpp             # Plain (direct elimination)
├── TransformationConstraintHandler.cpp # Transformation
└── AutoConstraintHandler.cpp    # Automatic selection
```

**femlabpy equivalent:** Lagrange and direct elimination
**Gap:** Penalty, transformation methods

### 3.5 Integration Algorithms

```
SRC/analysis/integrator/
├── StaticIntegrator.cpp              # Base static
├── IncrementalIntegrator.cpp         # Incremental
├── TransientIntegrator.cpp           # Base transient
├── ArcLength.cpp                     # Arc-length
├── ArcLength1.cpp                    # Arc-length variant 1
├── DisplacementControl.cpp           # Displacement control
├── LoadControl.cpp                   # Load control
├── MinUnbalDispNorm.cpp              # Min unbalanced disp
├── Newmark.cpp                       # Newmark-beta
├── Newmark1.cpp                      # Newmark variant 1
├── NewmarkExplicit.cpp               # Explicit Newmark
├── NewmarkHSFixedNumIter.cpp         # HS fixed iterations
├── NewmarkHSIncrLimit.cpp            # HS increment limit
├── NewmarkHSIncrReduct.cpp           # HS increment reduction
├── AlphaOS.cpp                       # Alpha-OS
├── AlphaOS_TP.cpp                    # Alpha-OS TP
├── AlphaOSGeneralized.cpp            # Generalized Alpha-OS
├── AlphaOSGeneralized_TP.cpp         # Generalized Alpha-OS TP
├── CentralDifference.cpp             # Central difference
├── CentralDifferenceAlternative.cpp  # Alternative CD
├── CentralDifferenceNoDamping.cpp    # CD without damping
├── Collocation.cpp                   # Collocation
├── CollocationHSFixedNumIter.cpp     # HS Collocation
├── CollocationHSIncrLimit.cpp        # HS increment limit
├── CollocationHSIncrReduct.cpp       # HS reduction
├── GeneralizedAlpha.cpp              # Generalized-alpha
├── HHT.cpp                           # HHT-alpha
├── HHT1.cpp                          # HHT variant 1
├── HHTExplicit.cpp                   # Explicit HHT
├── HHTExplicit_TP.cpp                # HHT Explicit TP
├── HHTGeneralized.cpp                # Generalized HHT
├── HHTGeneralized_TP.cpp             # Generalized HHT TP
├── HHTHSFixedNumIter.cpp             # HHT HS fixed
├── HHTHSFixedNumIter_TP.cpp          # HHT HS TP
├── HHTHSIncrLimit.cpp                # HHT HS limit
├── HHTHSIncrLimit_TP.cpp             # HHT HS limit TP
├── HHTHSIncrReduct.cpp               # HHT HS reduction
├── HHTHSIncrReduct_TP.cpp            # HHT HS reduction TP
├── TRBDF2.cpp                        # TR-BDF2
├── WilsonTheta.cpp                   # Wilson-theta
├── EQPath.cpp                        # Equilibrium path
└── HSConstraint.cpp                  # HS constraint
```

**femlabpy equivalent:** Newmark, Central Difference, HHT
**Gap:** Arc-length, displacement control, Wilson-theta, generalized-alpha, TRBDF2, collocation methods

---

## 4. CalculiX Source Code Analysis

### 4.1 Solver Architecture

CalculiX CCX is written primarily in Fortran with C interfaces. The solver supports:

- **Input Format:** ABAQUS-compatible input deck
- **Output Format:** frd (CalculiX native), VTK, dat
- **Solver Types:** SPOOLES, PARDISO, iterative (CG, BiCGSTAB)

### 4.2 Element Types

CalculiX supports the following element families:

#### 4.2.1 Continuum Elements (3D)
| Element | Nodes | Integration | Description |
|---------|-------|-------------|-------------|
| C3D4 | 4 | 1 point | Linear tetrahedron |
| C3D10 | 10 | 4 points | Quadratic tetrahedron |
| C3D10T | 10 | 4 points | Quadratic tet, thermal |
| C3D10I | 10 | 4 points | Improved quadratic tet |
| C3D8 | 8 | 8 points | Linear hexahedron |
| C3D8R | 8 | 1 point | Reduced integration hex |
| C3D8I | 8 | 8 points | Incompatible modes hex |
| C3D20 | 20 | 27 points | Quadratic hexahedron |
| C3D20R | 20 | 8 points | Reduced integration |
| C3D6 | 6 | 2 points | Linear wedge |
| C3D15 | 15 | 9 points | Quadratic wedge |

**femlabpy equivalent:** T4, H8
**Gap:** T10, H8R, H8I, H20, H20R, C3D6, C3D15

#### 4.2.2 Plane Elements (2D)
| Element | Nodes | Type | Description |
|---------|-------|------|-------------|
| CPS3 | 3 | PS | Linear triangle, plane stress |
| CPS4 | 4 | PS | Linear quad, plane stress |
| CPS4R | 4 | PS | Reduced integration |
| CPS6 | 6 | PS | Quadratic triangle |
| CPS8 | 8 | PS | Quadratic quad |
| CPS8R | 8 | PS | Reduced quadratic |
| CPE3 | 3 | PE | Linear triangle, plane strain |
| CPE4 | 4 | PE | Linear quad, plane strain |
| CPE4R | 4 | PE | Reduced integration |
| CPE6 | 6 | PE | Quadratic triangle |
| CPE8 | 8 | PE | Quadratic quad |
| CPE8R | 8 | PE | Reduced quadratic |
| CAX3 | 3 | AX | Axisymmetric triangle |
| CAX4 | 4 | AX | Axisymmetric quad |
| CAX6 | 6 | AX | Quadratic axisymmetric tri |
| CAX8 | 8 | AX | Quadratic axisymmetric quad |

**femlabpy equivalent:** T3, Q4 (plane stress/strain)
**Gap:** T6, Q8, reduced integration, axisymmetric

#### 4.2.3 Shell Elements
| Element | Nodes | Description |
|---------|-------|-------------|
| S3 | 3 | Linear triangle shell |
| S4 | 4 | Linear quad shell |
| S4R | 4 | Reduced integration |
| S6 | 6 | Quadratic triangle shell |
| S8 | 8 | Quadratic quad shell |
| S8R | 8 | Reduced quadratic |

**femlabpy equivalent:** NONE
**Gap:** Complete shell element family missing

#### 4.2.4 Beam Elements
| Element | Nodes | Description |
|---------|-------|-------------|
| B21 | 2 | 2D linear beam (Euler) |
| B22 | 3 | 2D quadratic beam |
| B31 | 2 | 3D linear beam (Euler) |
| B31R | 2 | 3D reduced integration |
| B32 | 3 | 3D quadratic beam |
| B32R | 3 | 3D reduced quadratic |

**femlabpy equivalent:** NONE
**Gap:** Complete beam element family missing

#### 4.2.5 Truss Elements
| Element | Nodes | Description |
|---------|-------|-------------|
| T2D2 | 2 | 2D truss |
| T3D2 | 2 | 3D truss |
| T3D3 | 3 | 3D quadratic truss |

**femlabpy equivalent:** Bar elements in `bars.py`
**Gap:** Quadratic truss

#### 4.2.6 Spring and Connector Elements
| Element | Description |
|---------|-------------|
| SPRING1 | Grounded spring |
| SPRING2 | Two-node spring |
| SPRINGA | Axial spring |
| DASHPOT1 | Grounded dashpot |
| DASHPOT2 | Two-node dashpot |
| DCOUP3D | Coupling element |
| GAP | Gap element |
| GAPUNI | Unidirectional gap |

**femlabpy equivalent:** NONE
**Gap:** Complete spring/connector family missing

#### 4.2.7 Mass and Rotary Inertia
| Element | Description |
|---------|-------------|
| MASS | Point mass |
| ROTARYI | Rotary inertia |

**femlabpy equivalent:** NONE
**Gap:** Point mass elements missing

### 4.3 Material Models

#### 4.3.1 Elastic Materials
- **ELASTIC**: Isotropic elastic
- **ELASTIC,TYPE=ORTHO**: Orthotropic elastic
- **ELASTIC,TYPE=ANISO**: Anisotropic elastic
- **ELASTIC,TYPE=ENGINEERING CONSTANTS**: Engineering notation

**femlabpy equivalent:** Isotropic only
**Gap:** Orthotropic, anisotropic materials

#### 4.3.2 Plastic Materials
- ***PLASTIC**: J2 plasticity
- ***PLASTIC,HARDENING=ISOTROPIC**: Isotropic hardening
- ***PLASTIC,HARDENING=KINEMATIC**: Kinematic hardening
- ***PLASTIC,HARDENING=COMBINED**: Combined hardening
- ***PLASTIC,HARDENING=JOHNSON COOK**: Johnson-Cook

**femlabpy equivalent:** Isotropic hardening von Mises
**Gap:** Kinematic hardening, combined hardening, Johnson-Cook

#### 4.3.3 Hyperelastic Materials
- ***HYPERELASTIC,MOONEY-RIVLIN**: Mooney-Rivlin
- ***HYPERELASTIC,POLYNOMIAL**: Polynomial form
- ***HYPERELASTIC,REDUCED POLYNOMIAL**: Reduced polynomial
- ***HYPERELASTIC,NEO HOOKE**: Neo-Hookean
- ***HYPERELASTIC,YEOH**: Yeoh model
- ***HYPERELASTIC,ARRUDA-BOYCE**: Arruda-Boyce
- ***HYPERELASTIC,OGDEN**: Ogden model

**femlabpy equivalent:** NONE
**Gap:** Complete hyperelastic library missing

#### 4.3.4 Viscoplastic Materials
- ***CREEP**: Creep model
- ***CREEP,LAW=NORTON**: Norton creep
- ***CREEP,LAW=USER**: User-defined creep
- ***RATE DEPENDENT**: Strain rate dependent

**femlabpy equivalent:** NONE
**Gap:** Creep and viscoplasticity missing

#### 4.3.5 Damage Models
- ***DAMAGE INITIATION**: Ductile, shear, FLD
- ***DAMAGE EVOLUTION**: Damage evolution

**femlabpy equivalent:** NONE
**Gap:** Damage mechanics missing

#### 4.3.6 User Materials
- ***USER MATERIAL**: UMAT interface

**femlabpy equivalent:** NONE
**Gap:** User material interface missing

### 4.4 Analysis Capabilities

#### 4.4.1 Static Analysis
- ***STATIC**: Linear/nonlinear static
- ***STATIC,NLGEOM**: Geometric nonlinearity
- ***STATIC,SOLVER=...**: Various solvers
- ***BUCKLE**: Linear buckling

**femlabpy equivalent:** Linear static
**Gap:** Geometric nonlinearity (general), buckling

#### 4.4.2 Dynamic Analysis
- ***DYNAMIC**: Implicit dynamics (Newmark)
- ***DYNAMIC,ALPHA=-0.05**: HHT-alpha
- ***MODAL DYNAMIC**: Modal superposition
- ***STEADY STATE DYNAMICS**: Frequency response
- ***RESPONSE SPECTRUM**: Response spectrum

**femlabpy equivalent:** Newmark, HHT, Central Difference
**Gap:** Modal superposition, steady-state, response spectrum

#### 4.4.3 Eigenvalue Analysis
- ***FREQUENCY**: Natural frequencies
- ***COMPLEX FREQUENCY**: Complex eigenvalues
- ***FREQUENCY,SOLVER=...**: Various solvers (ARPACK, SPOOLES)

**femlabpy equivalent:** Basic eigenvalue
**Gap:** Complex eigenvalues

#### 4.4.4 Thermal Analysis
- ***HEAT TRANSFER**: Steady/transient heat
- ***COUPLED TEMPERATURE-DISPLACEMENT**: Thermo-mechanical

**femlabpy equivalent:** NONE
**Gap:** Complete thermal analysis missing

#### 4.4.5 Coupled Analysis
- ***UNCOUPLED TEMPERATURE-DISPLACEMENT**: Sequentially coupled
- ***COUPLED TEMPERATURE-DISPLACEMENT**: Fully coupled

**femlabpy equivalent:** NONE
**Gap:** Coupled analysis missing

---

## 5. femlabpy Current State Analysis

### 5.1 Existing Elements

| Module | Elements | Coverage |
|--------|----------|----------|
| `elements/triangles.py` | T3 (CST) | Basic only |
| `elements/quads.py` | Q4 | Basic + elastoplastic |
| `elements/bars.py` | Truss | Basic + geometric NL |
| `elements/solids.py` | T4, H8 | Basic only |

**Total: 6 element types**

### 5.2 Existing Materials

| Module | Models | Coverage |
|--------|--------|----------|
| `materials/plasticity.py` | von Mises | Plane stress/strain |
| `materials/invariants.py` | Stress invariants | Helper functions |

**Total: 1 material model**

### 5.3 Existing Solvers

| Module | Solvers | Coverage |
|--------|---------|----------|
| `solvers.py` | `solve_nlbar` | Nonlinear truss |
| `solvers.py` | `solve_plastic` | Plasticity driver |
| `dynamics.py` | `solve_newmark` | Implicit Newmark |
| `dynamics.py` | `solve_hht` | HHT-alpha |
| `dynamics.py` | `solve_central_diff` | Explicit |
| `dynamics.py` | `solve_newmark_nl` | Nonlinear Newmark |
| `modal.py` | `solve_modal` | Eigenvalue |

**Total: 7 solver routines**

### 5.4 Existing I/O

| Module | Formats | Coverage |
|--------|---------|----------|
| `io/gmsh.py` | Gmsh 2.x/4.x | Mesh import |

**Total: 1 mesh format**

---

## 6. Gap Analysis Matrix

### 6.1 Element Gap Summary

| Category | OpenSees | CalculiX | femlabpy | Priority |
|----------|----------|----------|----------|----------|
| **1D Elements** | | | | |
| Truss | 4 types | 3 types | 1 type | Medium |
| Beam (2D) | 15+ types | 3 types | 0 | **HIGH** |
| Beam (3D) | 15+ types | 4 types | 0 | **HIGH** |
| Pipe | 2 types | N/A | 0 | Low |
| Cable | 2 types | N/A | 0 | Low |
| **2D Elements** | | | | |
| T3 | 1 type | 6 types | 1 type | Done |
| T6 | 1 type | 6 types | 0 | **HIGH** |
| Q4 | 9 types | 6 types | 1 type | Partial |
| Q8 | 1 type | 6 types | 0 | **HIGH** |
| Q9 | 1 type | N/A | 0 | Medium |
| Axisymmetric | N/A | 4 types | 0 | Medium |
| **3D Elements** | | | | |
| T4 | 1 type | 4 types | 1 type | Done |
| T10 | 2 types | 4 types | 0 | **HIGH** |
| H8 | 7 types | 4 types | 1 type | Partial |
| H20 | 3 types | 2 types | 0 | **HIGH** |
| H27 | 1 type | N/A | 0 | Medium |
| Wedge | N/A | 2 types | 0 | Low |
| **Shell Elements** | | | | |
| Triangular | 3 types | 2 types | 0 | **HIGH** |
| Quadrilateral | 10 types | 4 types | 0 | **HIGH** |
| **Special Elements** | | | | |
| Zero-length | 8 types | N/A | 0 | **HIGH** |
| Joint | 6 types | N/A | 0 | Medium |
| Bearing | 15 types | N/A | 0 | Medium |
| Contact | 5+ types | 2 types | 0 | Medium |
| Mass | N/A | 2 types | 0 | Medium |
| Spring | 1 type | 4 types | 0 | **HIGH** |

### 6.2 Material Gap Summary

| Category | OpenSees | CalculiX | femlabpy | Priority |
|----------|----------|----------|----------|----------|
| **Uniaxial Materials** | | | | |
| Steel | 17 types | 1 type | 0 | **HIGH** |
| Concrete | 26 types | 1 type | 0 | **HIGH** |
| Hysteretic | 24 types | N/A | 0 | **HIGH** |
| Damper | 11 types | 1 type | 0 | Medium |
| Gap/Contact | 7 types | 1 type | 0 | Medium |
| SMA | 4 types | N/A | 0 | Low |
| Isolator | 3 types | N/A | 0 | Low |
| Wrapper | 13 types | N/A | 0 | Medium |
| **nD Materials** | | | | |
| Elastic | 6 types | 4 types | 1 type | Partial |
| Plastic (J2) | 4 types | 4 types | 1 type | Partial |
| Drucker-Prager | 3 types | 1 type | 0 | **HIGH** |
| Mohr-Coulomb | N/A | 1 type | 0 | Medium |
| Soil (PDMY) | 6 types | N/A | 0 | Medium |
| Damage | 2 types | 1 type | 0 | Medium |
| **Hyperelastic** | | | | |
| Neo-Hookean | N/A | 1 type | 0 | **HIGH** |
| Mooney-Rivlin | N/A | 1 type | 0 | **HIGH** |
| Ogden | N/A | 1 type | 0 | Medium |
| **Section Models** | | | | |
| Elastic | 8 types | N/A | 0 | **HIGH** |
| Fiber | 12 types | N/A | 0 | **HIGH** |
| Layered | 3 types | N/A | 0 | Medium |

### 6.3 Analysis Gap Summary

| Category | OpenSees | CalculiX | femlabpy | Priority |
|----------|----------|----------|----------|----------|
| **Static Analysis** | | | | |
| Linear | ✅ | ✅ | ✅ | Done |
| Geometric NL | ✅ | ✅ | Partial | **HIGH** |
| Material NL | ✅ | ✅ | Partial | Partial |
| Buckling | N/A | ✅ | ❌ | **HIGH** |
| **Dynamic Analysis** | | | | |
| Newmark | ✅ | ✅ | ✅ | Done |
| HHT-alpha | ✅ | ✅ | ✅ | Done |
| Central Diff | ✅ | ✅ | ✅ | Done |
| Wilson-theta | ✅ | N/A | ❌ | Medium |
| Gen-alpha | ✅ | N/A | ❌ | Medium |
| Modal superposition | ✅ | ✅ | ❌ | **HIGH** |
| Response spectrum | ✅ | ✅ | ❌ | **HIGH** |
| Steady-state | ✅ | ✅ | ❌ | Medium |
| **Eigenvalue** | | | | |
| Standard | ✅ | ✅ | ✅ | Done |
| Complex | N/A | ✅ | ❌ | Low |
| **Solution Algorithms** | | | | |
| Newton-Raphson | ✅ | ✅ | ✅ | Done |
| Modified Newton | ✅ | ✅ | ❌ | Medium |
| BFGS | ✅ | N/A | ❌ | Medium |
| Line search | ✅ | ✅ | ❌ | **HIGH** |
| Arc-length | ✅ | N/A | ❌ | **HIGH** |
| **Thermal** | | | | |
| Steady-state | N/A | ✅ | ❌ | **HIGH** |
| Transient | N/A | ✅ | ❌ | **HIGH** |
| Coupled | N/A | ✅ | ❌ | Medium |

---

## 7. Implementation Roadmap

### Phase 1: Core Structural Elements (3-4 months)

**Objective:** Add essential beam and shell elements to enable structural frame analysis.

#### Phase 1.1: 2D Beam Elements (4-6 weeks)
```
Priority: CRITICAL
Reference: OpenSees SRC/element/beam2d/, SRC/element/elasticBeamColumn/
```

##### Task 1.1.1: Euler-Bernoulli Beam 2D
- [ ] Create `src/femlabpy/elements/beams2d.py`
- [ ] Implement `keEulerBernoulli2d()` - 6x6 stiffness matrix
- [ ] Implement `meEulerBernoulli2d()` - 6x6 consistent mass matrix
- [ ] Implement `meEulerBernoulli2d_lumped()` - Lumped mass
- [ ] Implement `qeEulerBernoulli2d()` - Internal forces
- [ ] Add coordinate transformation (local to global)
- [ ] Add distributed load equivalent nodal forces
- [ ] Write comprehensive tests with analytical solutions
- [ ] Benchmark against OpenSees elasticBeamColumn

##### Task 1.1.2: Timoshenko Beam 2D
- [ ] Implement `keTimoshenko2d()` - 6x6 with shear deformation
- [ ] Implement reduced integration for shear locking prevention
- [ ] Add shear correction factor calculation
- [ ] Write tests comparing thick vs thin beams
- [ ] Validate against CalculiX B21/B22

##### Task 1.1.3: Beam with P-Delta
- [ ] Implement geometric stiffness matrix `kg_beam2d()`
- [ ] Add iterative solution for geometric nonlinearity
- [ ] Write buckling examples
- [ ] Benchmark against OpenSees beam2d03

##### Task 1.1.4: Assembly Functions
- [ ] Implement `kbeam2d()` - Global assembly
- [ ] Implement `mbeam2d()` - Mass assembly
- [ ] Implement `qbeam2d()` - Force assembly
- [ ] Add beam load functions (distributed, point)

#### Phase 1.2: 3D Beam Elements (4-6 weeks)
```
Priority: HIGH
Reference: OpenSees SRC/element/beam3d/, SRC/element/elasticBeamColumn/
```

##### Task 1.2.1: Euler-Bernoulli Beam 3D
- [ ] Create `src/femlabpy/elements/beams3d.py`
- [ ] Implement `keEulerBernoulli3d()` - 12x12 stiffness matrix
- [ ] Implement `meEulerBernoulli3d()` - 12x12 mass matrix
- [ ] Add torsion stiffness (St. Venant)
- [ ] Add coordinate transformation (3D rotation matrix)
- [ ] Handle arbitrary cross-section orientation
- [ ] Write comprehensive tests

##### Task 1.2.2: Timoshenko Beam 3D
- [ ] Implement `keTimoshenko3d()` - 12x12 with shear
- [ ] Add 2x2 shear correction factor matrix
- [ ] Implement reduced integration
- [ ] Validate against CalculiX B31/B32

##### Task 1.2.3: Beam Orientation and Releases
- [ ] Implement `beam_orientation_vector()`
- [ ] Add member end releases (pins)
- [ ] Add rigid offsets
- [ ] Test with frame examples

#### Phase 1.3: Shell Elements (6-8 weeks)
```
Priority: HIGH
Reference: OpenSees SRC/element/shell/, CalculiX S4/S8
```

##### Task 1.3.1: Flat Shell Theory Review
- [ ] Document DKT (Discrete Kirchhoff Triangle) formulation
- [ ] Document MITC4 formulation
- [ ] Compare Reissner-Mindlin vs Kirchhoff theories
- [ ] Decide on implementation approach

##### Task 1.3.2: DKT Shell Triangle
- [ ] Create `src/femlabpy/elements/shells.py`
- [ ] Implement DKT bending stiffness (9 DOF reduced to 9 rotation DOF)
- [ ] Implement membrane stiffness (CST-based)
- [ ] Combine into full shell (18 DOF → reduced)
- [ ] Add drilling DOF stabilization
- [ ] Implement mass matrix
- [ ] Write patch tests

##### Task 1.3.3: MITC4 Shell Quad
- [ ] Implement MITC4 formulation
- [ ] Add assumed shear strain interpolation
- [ ] Implement 2x2 Gauss integration
- [ ] Add drilling DOF
- [ ] Write patch tests
- [ ] Benchmark against OpenSees ShellMITC4

##### Task 1.3.4: Shell Assembly and Post-processing
- [ ] Implement `kshell()` global assembly
- [ ] Add shell stress recovery (top/middle/bottom)
- [ ] Add shell resultants (Mx, My, Mxy, Vx, Vy, Nx, Ny, Nxy)
- [ ] Implement visualization helpers

#### Phase 1.4: Spring and Link Elements (2-3 weeks)
```
Priority: HIGH
Reference: OpenSees SRC/element/zeroLength/
```

##### Task 1.4.1: Zero-Length Spring
- [ ] Create `src/femlabpy/elements/springs.py`
- [ ] Implement `ke_spring()` - DOF-to-DOF spring
- [ ] Support translational and rotational DOFs
- [ ] Add multi-DOF springs
- [ ] Write basic tests

##### Task 1.4.2: Gap Element
- [ ] Implement compression-only gap
- [ ] Implement tension-only gap
- [ ] Add contact stiffness
- [ ] Test gap closure behavior

##### Task 1.4.3: Two-Node Link
- [ ] Implement general link element
- [ ] Support arbitrary transformation
- [ ] Add damper capability
- [ ] Test isolation systems

---

### Phase 2: Material Library Expansion (4-6 months)

**Objective:** Build comprehensive material library for structural analysis.

#### Phase 2.1: Steel Materials (4-6 weeks)
```
Priority: CRITICAL
Reference: OpenSees SRC/material/uniaxial/Steel01.cpp, Steel02.cpp
```

##### Task 2.1.1: Steel01 - Bilinear with Kinematic Hardening
- [ ] Create `src/femlabpy/materials/steel.py`
- [ ] Implement bilinear backbone curve
- [ ] Implement kinematic hardening (Bauschinger effect)
- [ ] Add isotropic hardening option
- [ ] Write cyclic loading tests
- [ ] Validate against OpenSees Steel01

##### Task 2.1.2: Steel02 - Giuffré-Menegotto-Pinto
- [ ] Implement GMP model equations
- [ ] Add curvature parameter R
- [ ] Implement hardening parameters (a1, a2, a3, a4)
- [ ] Add strain hardening
- [ ] Write monotonic and cyclic tests
- [ ] Validate against OpenSees Steel02

##### Task 2.1.3: Reinforcing Steel
- [ ] Implement ReinforcingSteel model
- [ ] Add bar buckling
- [ ] Add fatigue effects
- [ ] Test with concrete sections

##### Task 2.1.4: Steel with Fracture
- [ ] Implement ductile fracture model
- [ ] Add damage index tracking
- [ ] Add strength/stiffness degradation
- [ ] Test fracture initiation

#### Phase 2.2: Concrete Materials (6-8 weeks)
```
Priority: CRITICAL
Reference: OpenSees SRC/material/uniaxial/Concrete01.cpp, Concrete02.cpp
```

##### Task 2.2.1: Concrete01 - Kent-Scott-Park
- [ ] Create `src/femlabpy/materials/concrete.py`
- [ ] Implement parabolic compression curve
- [ ] Implement linear descending branch
- [ ] Add zero tensile strength
- [ ] Implement unloading/reloading rules
- [ ] Validate against OpenSees Concrete01

##### Task 2.2.2: Concrete02 - Linear Tension Softening
- [ ] Add linear tension softening branch
- [ ] Implement tension stiffening
- [ ] Add post-cracking behavior
- [ ] Write tension tests
- [ ] Validate against OpenSees Concrete02

##### Task 2.2.3: Confined Concrete
- [ ] Implement Mander model for confined concrete
- [ ] Add confinement effectiveness factor
- [ ] Support circular and rectangular sections
- [ ] Test with various confinement levels

##### Task 2.2.4: Concrete04 - Popovics
- [ ] Implement Popovics compression model
- [ ] Add Tsai tension model
- [ ] Implement cyclic rules
- [ ] Add crack closure

##### Task 2.2.5: Concrete with Damage
- [ ] Implement damage parameter tracking
- [ ] Add stiffness degradation
- [ ] Add strength degradation
- [ ] Test damage evolution

#### Phase 2.3: Hysteretic Materials (4-6 weeks)
```
Priority: HIGH
Reference: OpenSees SRC/material/uniaxial/HystereticMaterial.cpp
```

##### Task 2.3.1: Generic Hysteretic
- [ ] Create `src/femlabpy/materials/hysteretic.py`
- [ ] Implement multilinear backbone
- [ ] Add pinching parameters
- [ ] Add damage parameters (ductility, energy, unloading)
- [ ] Implement reloading rules
- [ ] Write cyclic tests

##### Task 2.3.2: Bouc-Wen Model
- [ ] Implement Bouc-Wen differential equation
- [ ] Add α, β, γ, n parameters
- [ ] Implement time integration for hysteretic variable
- [ ] Test under various loadings

##### Task 2.3.3: IMK Models
- [ ] Implement Modified Ibarra-Medina-Krawinkler model
- [ ] Add peak-oriented version
- [ ] Add pinching version
- [ ] Implement deterioration rules
- [ ] Validate for column hinges

##### Task 2.3.4: Pinching4
- [ ] Implement Pinching4 backbone
- [ ] Add four pinching parameters
- [ ] Add four damage parameters
- [ ] Test wall behavior

#### Phase 2.4: Fiber Sections (4-6 weeks)
```
Priority: HIGH
Reference: OpenSees SRC/material/section/FiberSection2d.cpp
```

##### Task 2.4.1: Fiber Section Framework
- [ ] Create `src/femlabpy/materials/sections.py`
- [ ] Define fiber data structure (y, z, area, material)
- [ ] Implement section force-deformation
- [ ] Add 2D and 3D versions
- [ ] Test with uniform fibers

##### Task 2.4.2: Section Discretization
- [ ] Implement rectangular patch discretization
- [ ] Implement circular/sector patch
- [ ] Implement straight/circular reinforcing layers
- [ ] Add I-section, T-section templates

##### Task 2.4.3: Fiber Section Integration
- [ ] Implement section state determination
- [ ] Add tangent stiffness computation
- [ ] Implement stress recovery
- [ ] Test moment-curvature response

##### Task 2.4.4: Section Aggregator
- [ ] Implement section aggregation
- [ ] Add shear flexibility
- [ ] Add torsion flexibility
- [ ] Test complete sections

#### Phase 2.5: Hyperelastic Materials (3-4 weeks)
```
Priority: MEDIUM
Reference: CalculiX hyperelastic models
```

##### Task 2.5.1: Neo-Hookean
- [ ] Create `src/femlabpy/materials/hyperelastic.py`
- [ ] Implement strain energy function
- [ ] Derive stress from energy
- [ ] Derive consistent tangent
- [ ] Test large deformation

##### Task 2.5.2: Mooney-Rivlin
- [ ] Implement 2-parameter Mooney-Rivlin
- [ ] Add incompressibility penalty
- [ ] Derive tangent modulus
- [ ] Test rubber materials

##### Task 2.5.3: Ogden Model
- [ ] Implement Ogden strain energy
- [ ] Support multiple terms
- [ ] Implement tangent
- [ ] Test against published data

#### Phase 2.6: Soil Materials (4-6 weeks)
```
Priority: MEDIUM
Reference: OpenSees SRC/material/nD/PressureDependMultiYield.cpp
```

##### Task 2.6.1: Drucker-Prager 3D
- [ ] Create `src/femlabpy/materials/geomaterials.py`
- [ ] Implement DP yield surface
- [ ] Add associative flow rule
- [ ] Implement return mapping
- [ ] Test triaxial compression

##### Task 2.6.2: Mohr-Coulomb
- [ ] Implement MC yield surface
- [ ] Handle corner regions
- [ ] Add dilation angle
- [ ] Test slope stability

##### Task 2.6.3: PDMY (Simplified)
- [ ] Implement simplified pressure-dependent model
- [ ] Add multi-surface plasticity (reduced)
- [ ] Test cyclic behavior
- [ ] Compare with OpenSees

---

### Phase 3: Advanced Analysis (3-4 months)

**Objective:** Expand analysis capabilities for production use.

#### Phase 3.1: Geometric Nonlinearity (4-6 weeks)
```
Priority: CRITICAL
Reference: OpenSees corotational formulation
```

##### Task 3.1.1: Corotational Framework
- [ ] Create `src/femlabpy/analysis/corotational.py`
- [ ] Implement corotational coordinate system
- [ ] Derive transformation matrices
- [ ] Implement tangent stiffness transformation
- [ ] Test rigid body modes

##### Task 3.1.2: Corotational Truss
- [ ] Add corotational truss element
- [ ] Test snap-through behavior
- [ ] Validate against current nlbar
- [ ] Compare performance

##### Task 3.1.3: Corotational Beam 2D
- [ ] Implement corotational beam transformation
- [ ] Handle finite rotations
- [ ] Test frame buckling
- [ ] Validate against OpenSees

##### Task 3.1.4: Corotational Beam 3D
- [ ] Extend to 3D rotations (quaternions or rotation vectors)
- [ ] Handle 3D rigid body modes
- [ ] Test space frame buckling
- [ ] Benchmark results

##### Task 3.1.5: Updated Lagrangian for Solids
- [ ] Implement total Lagrangian formulation
- [ ] Add geometric stiffness for solids
- [ ] Test large deformation problems
- [ ] Validate against CalculiX NLGEOM

#### Phase 3.2: Arc-Length Methods (3-4 weeks)
```
Priority: HIGH
Reference: OpenSees SRC/analysis/integrator/ArcLength.cpp
```

##### Task 3.2.1: Spherical Arc-Length
- [ ] Create `src/femlabpy/analysis/arclength.py`
- [ ] Implement spherical constraint
- [ ] Add predictor step
- [ ] Implement corrector iterations
- [ ] Test snap-through

##### Task 3.2.2: Cylindrical Arc-Length
- [ ] Implement cylindrical constraint
- [ ] Add displacement-based constraint
- [ ] Test snap-back behavior

##### Task 3.2.3: Arc-Length with Line Search
- [ ] Add line search to arc-length
- [ ] Improve convergence
- [ ] Test difficult problems

##### Task 3.2.4: Displacement Control
- [ ] Implement displacement control loading
- [ ] Add multi-point displacement control
- [ ] Test with material softening

#### Phase 3.3: Buckling Analysis (3-4 weeks)
```
Priority: HIGH
Reference: CalculiX *BUCKLE
```

##### Task 3.3.1: Linear Buckling
- [ ] Create `src/femlabpy/analysis/buckling.py`
- [ ] Implement linearized buckling eigenvalue
- [ ] Form geometric stiffness at reference load
- [ ] Solve generalized eigenvalue problem
- [ ] Extract critical load factors

##### Task 3.3.2: Buckling Mode Shapes
- [ ] Post-process eigenvectors
- [ ] Visualize buckling modes
- [ ] Add mode normalization

##### Task 3.3.3: Post-Buckling Analysis
- [ ] Use buckling mode as imperfection
- [ ] Run nonlinear analysis
- [ ] Track equilibrium path
- [ ] Compare linear vs nonlinear buckling

#### Phase 3.4: Response Spectrum Analysis (3-4 weeks)
```
Priority: HIGH
Reference: OpenSees, CalculiX *RESPONSE SPECTRUM
```

##### Task 3.4.1: Response Spectrum Framework
- [ ] Create `src/femlabpy/analysis/spectrum.py`
- [ ] Implement spectrum data structure
- [ ] Add standard spectrum shapes (EC8, ASCE 7)
- [ ] Support user-defined spectra

##### Task 3.4.2: Modal Combination
- [ ] Implement SRSS combination
- [ ] Implement CQC combination
- [ ] Implement ABSSUM (for reference)
- [ ] Compare combination methods

##### Task 3.4.3: Multi-Direction Spectrum
- [ ] Support 3-direction input
- [ ] Implement SRSS for directions
- [ ] Implement 30% rule (100-30-30)
- [ ] Test 3D structures

##### Task 3.4.4: Response Extraction
- [ ] Extract peak displacements
- [ ] Extract peak forces
- [ ] Extract base shear
- [ ] Generate design envelopes

#### Phase 3.5: Modal Superposition (2-3 weeks)
```
Priority: MEDIUM
Reference: CalculiX *MODAL DYNAMIC
```

##### Task 3.5.1: Modal Transformation
- [ ] Create `src/femlabpy/analysis/modal_superposition.py`
- [ ] Transform to modal coordinates
- [ ] Solve decoupled SDOF equations
- [ ] Transform back to physical

##### Task 3.5.2: Modal Damping
- [ ] Support modal damping ratios per mode
- [ ] Implement Rayleigh damping in modal space
- [ ] Compare with direct integration

##### Task 3.5.3: Mode Truncation
- [ ] Implement mode truncation criteria
- [ ] Add static correction (residual modes)
- [ ] Validate against full solution

#### Phase 3.6: Thermal Analysis (4-6 weeks)
```
Priority: MEDIUM
Reference: CalculiX *HEAT TRANSFER
```

##### Task 3.6.1: Steady-State Heat Transfer
- [ ] Create `src/femlabpy/analysis/thermal.py`
- [ ] Implement conductivity matrix assembly
- [ ] Add temperature boundary conditions
- [ ] Add heat flux boundary conditions
- [ ] Test 2D and 3D problems

##### Task 3.6.2: Transient Heat Transfer
- [ ] Implement heat capacity matrix
- [ ] Add time integration (backward Euler)
- [ ] Test transient problems
- [ ] Compare with analytical solutions

##### Task 3.6.3: Thermal Elements
- [ ] Add thermal T3/T6 elements
- [ ] Add thermal Q4/Q8 elements
- [ ] Add thermal T4/T10 elements
- [ ] Add thermal H8/H20 elements

##### Task 3.6.4: Thermo-Mechanical Coupling
- [ ] Implement sequentially coupled analysis
- [ ] Add thermal strain computation
- [ ] Test thermal stress problems
- [ ] Validate against CalculiX

---

### Phase 4: High-Performance Computing (2-3 months)

**Objective:** Improve performance for large-scale problems.

#### Phase 4.1: Sparse Matrix Assembly (2-3 weeks)
```
Priority: HIGH
```

##### Task 4.1.1: COO Assembly
- [ ] Create `src/femlabpy/assembly/sparse.py`
- [ ] Implement COO triplet assembly
- [ ] Add duplicate handling
- [ ] Convert to CSR/CSC

##### Task 4.1.2: Direct CSR Assembly
- [ ] Implement symbolic assembly (structure)
- [ ] Implement numeric assembly (values)
- [ ] Benchmark against COO

##### Task 4.1.3: Assembly Optimization
- [ ] Use numba JIT compilation
- [ ] Add parallel assembly
- [ ] Profile and optimize

#### Phase 4.2: Sparse Solvers (2-3 weeks)
```
Priority: HIGH
```

##### Task 4.2.1: Iterative Solvers
- [ ] Add CG with preconditioner
- [ ] Add GMRES for non-symmetric
- [ ] Implement diagonal preconditioning
- [ ] Implement ILU preconditioning

##### Task 4.2.2: Direct Solver Interfaces
- [ ] Interface with CHOLMOD (SuiteSparse)
- [ ] Interface with UMFPACK
- [ ] Add solver selection logic

##### Task 4.2.3: Solver Benchmarking
- [ ] Create benchmark suite
- [ ] Compare solver performance
- [ ] Document best practices

#### Phase 4.3: Parallel Computing (4-6 weeks)
```
Priority: MEDIUM
```

##### Task 4.3.1: Shared Memory Parallelism
- [ ] Add numba parallel loops
- [ ] Add ThreadPoolExecutor usage
- [ ] Profile parallel efficiency

##### Task 4.3.2: GPU Acceleration (Optional)
- [ ] Investigate CuPy integration
- [ ] Implement GPU assembly
- [ ] Implement GPU solver
- [ ] Benchmark GPU vs CPU

##### Task 4.3.3: Distributed Memory (Optional)
- [ ] Investigate mpi4py integration
- [ ] Implement domain decomposition
- [ ] Add distributed solver

---

### Phase 5: Advanced Elements (3-4 months)

**Objective:** Complete element library for specialized applications.

#### Phase 5.1: Higher-Order Elements (4-6 weeks)

##### Task 5.1.1: T6 Triangle
- [ ] Implement 6-node triangle shape functions
- [ ] Add 3-point Gauss integration
- [ ] Implement stress recovery
- [ ] Patch tests

##### Task 5.1.2: Q8 Quadrilateral
- [ ] Implement 8-node serendipity shape functions
- [ ] Add 3x3 Gauss integration
- [ ] Implement reduced integration (2x2)
- [ ] Test hourglass control

##### Task 5.1.3: T10 Tetrahedron
- [ ] Implement 10-node tet shape functions
- [ ] Add 4-point integration
- [ ] Test locking behavior
- [ ] Compare with linear tet

##### Task 5.1.4: H20 Hexahedron
- [ ] Implement 20-node serendipity shape functions
- [ ] Add 3x3x3 full integration
- [ ] Add 2x2x2 reduced integration
- [ ] Test bending-dominated problems

#### Phase 5.2: Specialized Elements (4-6 weeks)

##### Task 5.2.1: Axisymmetric Elements
- [ ] Create `src/femlabpy/elements/axisymmetric.py`
- [ ] Implement axisymmetric T3
- [ ] Implement axisymmetric Q4
- [ ] Add body force integration
- [ ] Test pressure vessel

##### Task 5.2.2: Interface/Cohesive Elements
- [ ] Implement 4-node interface element
- [ ] Add traction-separation law
- [ ] Implement damage mechanics
- [ ] Test delamination

##### Task 5.2.3: Enhanced Strain Elements
- [ ] Implement Q1E4 (Q4 with enhanced strain)
- [ ] Add incompatible mode shape functions
- [ ] Test bending behavior
- [ ] Compare with standard Q4

##### Task 5.2.4: Reduced Integration Elements
- [ ] Implement Q4R with hourglass control
- [ ] Implement H8R with hourglass control
- [ ] Add Flanagan-Belytschko stabilization
- [ ] Test explicit dynamics

#### Phase 5.3: Contact Elements (4-6 weeks)

##### Task 5.3.1: Node-to-Segment Contact
- [ ] Create `src/femlabpy/elements/contact.py`
- [ ] Implement gap detection
- [ ] Add penalty formulation
- [ ] Test Hertz contact

##### Task 5.3.2: Segment-to-Segment Contact
- [ ] Implement mortar method
- [ ] Add Lagrange multiplier formulation
- [ ] Test sliding contact

##### Task 5.3.3: Friction
- [ ] Implement Coulomb friction
- [ ] Add regularized friction
- [ ] Test stick-slip behavior

---

## 8. Detailed Implementation Tasks

### 8.1 Element Implementation Template

For each new element, follow this standard workflow:

```markdown
#### Element: [NAME]

**Reference Implementation:** [OpenSees/CalculiX file]

**Mathematical Formulation:**
- Shape functions: [describe]
- Integration scheme: [Gauss points]
- DOFs per node: [number]
- Total DOFs: [number]

**Implementation Checklist:**
- [ ] Create element module file
- [ ] Implement shape functions `N_[element](xi, eta, [zeta])`
- [ ] Implement shape function derivatives `dN_[element](xi, eta, [zeta])`
- [ ] Implement Jacobian computation
- [ ] Implement B-matrix assembly
- [ ] Implement stiffness matrix `ke[element]()`
- [ ] Implement mass matrix `me[element]()`
- [ ] Implement internal force vector `qe[element]()`
- [ ] Implement stress recovery
- [ ] Implement global assembly function `k[element]()`
- [ ] Implement global mass assembly `m[element]()`
- [ ] Implement global force assembly `q[element]()`
- [ ] Write unit tests (shape function values)
- [ ] Write patch tests (constant stress)
- [ ] Write benchmark tests (analytical solutions)
- [ ] Write documentation (docstrings)
- [ ] Update API reference
- [ ] Add usage examples

**Validation Tests:**
- [ ] Patch test: constant stress field
- [ ] Eigenvalue test: rigid body modes
- [ ] Benchmark: [specific analytical problem]
- [ ] Comparison: OpenSees/CalculiX results

**Files to Create/Modify:**
- `src/femlabpy/elements/[module].py` (new)
- `src/femlabpy/elements/__init__.py` (update)
- `tests/test_[module].py` (new)
- `docs/api.rst` (update)
- `docs/manual/ch[X]_[topic].md` (update)
```

### 8.2 Material Implementation Template

```markdown
#### Material: [NAME]

**Reference Implementation:** [OpenSees/CalculiX file]

**Mathematical Formulation:**
- Yield criterion: [equation]
- Flow rule: [associative/non-associative]
- Hardening: [isotropic/kinematic/combined]
- State variables: [list]

**Implementation Checklist:**
- [ ] Create material module file
- [ ] Implement yield function `f([material])`
- [ ] Implement yield function gradient `df_dsigma([material])`
- [ ] Implement plastic flow direction `n([material])`
- [ ] Implement hardening modulus `H([material])`
- [ ] Implement return mapping `stress_[material]()`
- [ ] Implement consistent tangent `tangent_[material]()`
- [ ] Implement state variable initialization
- [ ] Implement state variable update
- [ ] Write unit tests (elastic response)
- [ ] Write unit tests (plastic response)
- [ ] Write cyclic loading tests
- [ ] Write integration accuracy tests
- [ ] Write documentation
- [ ] Add examples

**Validation Tests:**
- [ ] Uniaxial tension
- [ ] Uniaxial compression
- [ ] Biaxial loading
- [ ] Cyclic loading
- [ ] Comparison with OpenSees/CalculiX
```

### 8.3 Analysis Implementation Template

```markdown
#### Analysis Type: [NAME]

**Reference Implementation:** [OpenSees/CalculiX file]

**Mathematical Formulation:**
- Governing equation: [equation]
- Discretization: [spatial/temporal]
- Solution strategy: [Newton/arc-length/etc.]

**Implementation Checklist:**
- [ ] Create analysis module file
- [ ] Implement predictor step
- [ ] Implement corrector iterations
- [ ] Implement convergence check
- [ ] Implement solution update
- [ ] Implement result storage
- [ ] Write tests with known solutions
- [ ] Write documentation
- [ ] Add examples

**Validation Tests:**
- [ ] Linear problem (exact solution)
- [ ] Nonlinear problem (published results)
- [ ] Comparison with OpenSees/CalculiX
```

---

## 9. Testing and Validation Strategy

### 9.1 Test Categories

#### 9.1.1 Unit Tests
- Shape function values at specific points
- Gauss point locations and weights
- Matrix symmetry and positive definiteness
- Individual material response

#### 9.1.2 Patch Tests
- Constant stress field recovery
- Rigid body mode preservation
- Equilibrium satisfaction

#### 9.1.3 Benchmark Tests
- Analytical solutions (beam deflection, plate bending)
- Published results (NAFEMS benchmarks)
- Textbook examples

#### 9.1.4 Comparison Tests
- OpenSees: Import TCL model, compare results
- CalculiX: Import INP file, compare results
- Tolerance: 1% for displacements, 5% for stresses

### 9.2 Benchmark Problems

#### Structural Benchmarks
1. Cantilever beam tip deflection
2. Simply supported beam natural frequencies
3. Circular plate center deflection
4. Scordelis-Lo roof (shell benchmark)
5. Cook's membrane (incompressibility)
6. MacNeal-Harder beam (shear locking)
7. Euler column buckling
8. Williams' eigenvalue problem (singularity)

#### Dynamic Benchmarks
1. SDOF spring-mass (Newmark accuracy)
2. MDOF frame (mode shapes)
3. Seismic building (time history)
4. Impact problem (central difference)
5. Earthquake engineering examples

#### Nonlinear Benchmarks
1. Two-bar snap-through truss
2. Williams toggle frame
3. Lee frame (geometric NL)
4. Elastoplastic plate (J2 plasticity)
5. Concrete beam (material NL)

### 9.3 Continuous Integration

```yaml
# .github/workflows/tests.yml additions
jobs:
  benchmark:
    runs-on: ubuntu-latest
    steps:
      - name: Run OpenSees comparison
        run: python -m pytest tests/benchmarks/opensees/ -v
      
      - name: Run CalculiX comparison
        run: python -m pytest tests/benchmarks/calculix/ -v
      
      - name: Generate comparison report
        run: python scripts/generate_benchmark_report.py
```

---

## 10. References

### 10.1 OpenSees Documentation
- [OpenSees Wiki](https://opensees.berkeley.edu/wiki/)
- [OpenSees GitHub](https://github.com/OpenSees/OpenSees)
- [OpenSeesPy Documentation](https://openseespydoc.readthedocs.io/)
- McKenna, F. (2011). "OpenSees: A Framework for Earthquake Engineering Simulation." *Computing in Science & Engineering*.

### 10.2 CalculiX Documentation
- [CalculiX Website](http://www.calculix.de/)
- [CalculiX CCX Manual](http://www.dhondt.de/ccx_2.21.pdf)
- Dhondt, G. (2004). *The Finite Element Method for Three-Dimensional Thermomechanical Applications*. Wiley.

### 10.3 Finite Element Theory
- Bathe, K.J. (2014). *Finite Element Procedures*. 2nd Edition.
- Hughes, T.J.R. (2000). *The Finite Element Method*. Dover.
- Zienkiewicz, Taylor, Zhu (2013). *The Finite Element Method*. 7th Edition.
- Cook et al. (2002). *Concepts and Applications of Finite Element Analysis*. 4th Edition.

### 10.4 Structural Dynamics
- Chopra, A.K. (2017). *Dynamics of Structures*. 5th Edition.
- Clough, R.W., Penzien, J. (2003). *Dynamics of Structures*. 3rd Edition.
- Newmark, N.M. (1959). "A Method of Computation for Structural Dynamics." *ASCE*.

### 10.5 Nonlinear Analysis
- Crisfield, M.A. (1991, 1997). *Non-linear Finite Element Analysis of Solids and Structures*. Vols 1 & 2.
- Belytschko, Liu, Moran (2014). *Nonlinear Finite Elements for Continua and Structures*. 2nd Edition.
- de Souza Neto, Perić, Owen (2008). *Computational Methods for Plasticity*. Wiley.

### 10.6 Material Modeling
- Lemaitre, J., Chaboche, J.L. (1990). *Mechanics of Solid Materials*. Cambridge.
- Simo, J.C., Hughes, T.J.R. (1998). *Computational Inelasticity*. Springer.
- Chen, W.F. (1982). *Plasticity in Reinforced Concrete*. McGraw-Hill.

### 10.7 OpenSees Material References
- Steel02: Filippou, F.C. et al. (1983). "Effects of Bond Deterioration on Hysteretic Behavior of Reinforced Concrete Joints."
- Concrete02: Scott, B.D. et al. (1982). "Stress-Strain Behavior of Concrete Confined by Overlapping Hoops at Low and High Strain Rates."
- IMK Models: Ibarra, L.F. et al. (2005). "Hysteretic Models that Incorporate Strength and Stiffness Deterioration."
- PDMY: Yang, Z. et al. (2003). "Computational Model for Cyclic Mobility and Associated Shear Deformation."

---

## Appendix A: OpenSees Element Directory Listing

```
SRC/element/
├── CEqElement/
├── HUelements/
├── IGA/
├── LHMYS/
├── PFEMElement/
├── PML/
├── RockingBC/
├── UP-ucsd/
├── UWelements/
├── XMUelements/
├── absorbentBoundaries/
├── adapter/
├── beam2d/
├── beam3d/
├── beamWithHinges/
├── brick/
├── catenaryCable/
├── componentElement/
├── dispBeamColumn/
├── dispBeamColumnInt/
├── dmglib/
├── elasticBeamColumn/
├── elastomericBearing/
├── feap/
├── fmkPlanarTruss/
├── forceBeamColumn/
├── fourNodeQuad/
├── frictionBearing/
├── generic/
├── gradientInelasticBeamColumn/
├── joint/
├── masonry/
├── mefi/
├── mixedBeamColumn/
├── mvlem/
├── nonlinearBeamColumn/
├── pipe/
├── pyMacro/
├── shell/
├── surfaceLoad/
├── tetrahedron/
├── triangle/
├── truss/
├── twoNodeLink/
├── updatedLagrangianBeamColumn/
└── zeroLength/
```

---

## Appendix B: OpenSees Uniaxial Material Listing (170+ Materials)

```
SRC/material/uniaxial/
├── APDFMD.cpp
├── APDMD.cpp
├── APDVFD.cpp
├── ASDConcrete1DMaterial.cpp
├── ASDSteel1DMaterial.cpp
├── ASD_SMA_3K.cpp
├── AxialSp.cpp
├── AxialSpHD.cpp
├── BWBN.cpp
├── BackboneMaterial.cpp
├── BarSlipMaterial.cpp
├── Bilin.cpp
├── Bilin02.cpp
├── BilinearOilDamper.cpp
├── Bond_SP01.cpp
├── BoucWenInfill.cpp
├── BoucWenMaterial.cpp
├── BoucWenOriginal.cpp
├── BraceMaterial.cpp
├── CFSSSWP.cpp
├── CFSWSWP.cpp
├── CableMaterial.cpp
├── Cast.cpp
├── Concrete01.cpp
├── Concrete01WithSITC.cpp
├── Concrete02.cpp
├── Concrete02IS.cpp
├── Concrete02Thermal.cpp
├── Concrete04.cpp
├── Concrete05.cpp
├── Concrete06.cpp
├── Concrete07.cpp
├── ConcreteCM.cpp
├── ConcreteD.cpp
├── ConcreteECThermal.cpp
├── ConcreteSakaiKawashima.cpp
├── ConcreteZBH_fitted.cpp
├── ConcreteZBH_original.cpp
├── ConcreteZBH_smoothed.cpp
├── ConcretewBeta.cpp
├── ConfinedConcrete01.cpp
├── ContinuumUniaxial.cpp
├── CoulombDamperMaterial.cpp
├── CreepMaterial.cpp
├── CreepShrinkageACI209.cpp
├── CubicSpline.cpp
├── DamperMaterial.cpp
├── DegradingPinchedBW.cpp
├── DoddRestr.cpp
├── Dodd_Restrepo.cpp
├── DowelType.cpp
├── DrainMaterial.cpp
├── DuctileFracture.cpp
├── ECC01.cpp
├── ENTMaterial.cpp
├── EPPGapMaterial.cpp
├── Elastic2Material.cpp
├── ElasticBDMaterial.cpp
├── ElasticBilin.cpp
├── ElasticMaterial.cpp
├── ElasticMaterialThermal.cpp
├── ElasticMultiLinear.cpp
├── ElasticPPMaterial.cpp
├── ElasticPowerFunc.cpp
├── ExternalUniaxialMaterial.cpp
├── FRCC.cpp
├── FRPConfinedConcrete.cpp
├── FRPConfinedConcrete02.cpp
├── FatigueMaterial.cpp
├── FedeasMaterial.cpp
├── FlagShapeMaterial.cpp
├── GMG_CyclicReinforcedConcrete.cpp
├── GNGMaterial.cpp
├── HardeningMaterial.cpp
├── HardeningMaterial2.cpp
├── Hertzdamp.cpp
├── HookGap.cpp
├── HyperbolicGapMaterial.cpp
├── HystereticAsym.cpp
├── HystereticMaterial.cpp
├── HystereticPoly.cpp
├── HystereticSMMaterial.cpp
├── HystereticSmooth.cpp
├── IMKBilin.cpp
├── IMKPeakOriented.cpp
├── IMKPinching.cpp
├── ImpactMaterial.cpp
├── InitStrainMaterial.cpp
├── InitStressMaterial.cpp
├── JankowskiImpact.cpp
├── KikuchiAikenHDR.cpp
├── KikuchiAikenLRB.cpp
├── Masonry.cpp
├── Masonry2.cpp
├── Masonryt.cpp
├── MaterialState.cpp
├── Maxwell.cpp
├── MinMaxMaterial.cpp
├── ModIMKPeakOriented.cpp
├── ModIMKPeakOriented02.cpp
├── ModIMKPinching.cpp
├── ModIMKPinching02.cpp
├── MultiLinear.cpp
├── MultiplierMaterial.cpp
├── Neoprene.cpp
├── OOHystereticMaterial.cpp
├── OriginCentered.cpp
├── ParallelMaterial.cpp
├── PathIndependentMaterial.cpp
├── PenaltyMaterial.cpp
├── Pinching4Material.cpp
├── PipeMaterial.cpp
├── PrestressedSteelMaterial.cpp
├── QbSandCPT.cpp
├── RambergOsgoodSteel.cpp
├── Ratchet.cpp
├── ReinforcingSteel.cpp
├── ResilienceLow.cpp
├── ResilienceMaterialHR.cpp
├── SAWSMaterial.cpp
├── SLModel.cpp
├── SMAMaterial.cpp
├── SPSW02.cpp
├── SecantConcrete.cpp
├── SelfCenteringMaterial.cpp
├── SeriesMaterial.cpp
├── ShearPanelMaterial.cpp
├── SimpleFractureMaterial.cpp
├── SmoothPSConcrete.cpp
├── StainlessECThermal.cpp
├── Steel01.cpp
├── Steel01Thermal.cpp
├── Steel02.cpp
├── Steel02Fatigue.cpp
├── Steel02Thermal.cpp
├── Steel03.cpp
├── Steel2.cpp
├── Steel4.cpp
├── SteelBRB.cpp
├── SteelDRC.cpp
├── SteelECThermal.cpp
├── SteelFractureDI.cpp
├── SteelMP.cpp
├── SteelMPF.cpp
├── TDConcrete.cpp
├── TDConcreteEXP.cpp
├── TDConcreteMC10.cpp
├── TDConcreteMC10NL.cpp
├── TDConcreteNL.cpp
├── TensionOnlyMaterial.cpp
├── TriMatrix.cpp
├── Trilinwp.cpp
├── Trilinwp2.cpp
├── Trilinwpd.cpp
├── TzSandCPT.cpp
├── UVCuniaxial.cpp
├── UniaxialJ2Plasticity.cpp
├── UniaxialMaterial.cpp
├── ViscoelasticGap.cpp
├── ViscousDamper.cpp
├── ViscousMaterial.cpp
└── pyUCLA.cpp
```

---

## Appendix C: CalculiX Element Types Summary

| Category | Elements |
|----------|----------|
| **3D Continuum** | C3D4, C3D10, C3D10T, C3D8, C3D8R, C3D8I, C3D20, C3D20R, C3D6, C3D15 |
| **2D Plane Stress** | CPS3, CPS4, CPS4R, CPS6, CPS8, CPS8R |
| **2D Plane Strain** | CPE3, CPE4, CPE4R, CPE6, CPE8, CPE8R |
| **Axisymmetric** | CAX3, CAX4, CAX6, CAX8 |
| **Shell** | S3, S4, S4R, S6, S8, S8R |
| **Beam** | B21, B22, B31, B31R, B32, B32R |
| **Truss** | T2D2, T3D2, T3D3 |
| **Connector** | SPRING1, SPRING2, SPRINGA, DASHPOT1, DASHPOT2, DCOUP3D, GAP, GAPUNI |
| **Mass** | MASS, ROTARYI |

---

## Appendix D: femlabpy API Expansion Summary

### Current API (v0.6.0)
```python
# Elements
femlabpy.ket3e, kt3e, qet3e, qt3e       # T3 triangle
femlabpy.keq4e, kq4e, qeq4e, qq4e       # Q4 quad
femlabpy.kebar, kbar, qebar, qbar       # Truss
femlabpy.keT4e, kT4e, qeT4e, qT4e       # T4 tet
femlabpy.keh8e, kh8e, qeh8e, qh8e       # H8 hex

# Solvers
femlabpy.elastic()                       # Linear elastic
femlabpy.plastps(), plastpe()           # Plasticity
femlabpy.solve_newmark()                # Newmark
femlabpy.solve_hht()                    # HHT-alpha
femlabpy.solve_central_diff()           # Central difference
femlabpy.solve_modal()                  # Eigenvalue
```

### Proposed API Expansion (v1.0.0+)
```python
# New Elements
femlabpy.elements.beam2d.*              # 2D beams
femlabpy.elements.beam3d.*              # 3D beams
femlabpy.elements.shell.*               # Shells
femlabpy.elements.springs.*             # Springs/links
femlabpy.elements.T6, Q8, T10, H20.*   # Higher-order

# New Materials
femlabpy.materials.steel.*              # Steel models
femlabpy.materials.concrete.*           # Concrete models
femlabpy.materials.hysteretic.*         # Hysteretic models
femlabpy.materials.hyperelastic.*       # Rubber models
femlabpy.materials.sections.*           # Fiber sections

# New Analysis
femlabpy.analysis.arclength.*           # Arc-length
femlabpy.analysis.buckling.*            # Buckling
femlabpy.analysis.spectrum.*            # Response spectrum
femlabpy.analysis.thermal.*             # Thermal
femlabpy.analysis.corotational.*        # Geometric NL
```

---

## Document History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-04-05 | Copilot | Initial comprehensive analysis |

---

## Appendix E: Detailed Element Mathematical Formulations

### E.1 Euler-Bernoulli Beam Element (2D)

#### E.1.1 Kinematics

The Euler-Bernoulli beam assumes:
- Plane sections remain plane and perpendicular to the neutral axis
- Shear deformation is neglected
- Small strain theory applies

Displacement field:
```
u(x, y) = u₀(x) - y × dv/dx
v(x, y) = v(x)
```

Where:
- u₀(x) = axial displacement of neutral axis
- v(x) = transverse displacement
- y = distance from neutral axis

#### E.1.2 Shape Functions

For a 2-node beam element with DOFs [u₁, v₁, θ₁, u₂, v₂, θ₂]:

Axial shape functions (linear):
```
N₁ᵃ(ξ) = (1 - ξ) / 2
N₂ᵃ(ξ) = (1 + ξ) / 2
```

Transverse shape functions (Hermite cubics):
```
H₁(ξ) = (1 - 3ξ² + 2ξ³) / 4
H₂(ξ) = L × (ξ - 2ξ² + ξ³) / 8
H₃(ξ) = (1 + 3ξ² - 2ξ³) / 4
H₄(ξ) = L × (-ξ - 2ξ² + ξ³) / 8
```

Where:
- ξ = 2x/L - 1 (natural coordinate, -1 ≤ ξ ≤ 1)
- L = element length

#### E.1.3 Stiffness Matrix

Element stiffness matrix (6×6):
```
K_e = [K_aa  0    0   | K_ab  0    0  ]
      [0    K_bb K_bc | 0    K_bd K_be]
      [0    K_cb K_cc | 0    K_cd K_ce]
      [-----|--------|-----|-----------|
      [K_ba  0    0   | K_dd  0    0  ]
      [0    K_db K_dc | 0    K_ee K_ef]
      [0    K_eb K_ec | 0    K_fe K_ff]
```

Where (explicit formulas):
```
K_aa = EA/L,  K_ab = -EA/L
K_bb = 12EI/L³,  K_bc = 6EI/L²,  K_bd = -12EI/L³,  K_be = 6EI/L²
K_cc = 4EI/L,  K_cd = -6EI/L²,  K_ce = 2EI/L
K_dd = EA/L
K_ee = 12EI/L³,  K_ef = -6EI/L²
K_ff = 4EI/L
```

#### E.1.4 Mass Matrix

Consistent mass matrix:
```
M_e = ρAL/420 × [140   0     0    | 70    0     0   ]
               [0    156   22L   | 0     54   -13L ]
               [0    22L   4L²   | 0     13L  -3L² ]
               [-----|-----------|-----|------------|
               [70    0     0    | 140   0     0   ]
               [0     54   13L   | 0    156  -22L ]
               [0    -13L  -3L²  | 0   -22L   4L² ]
```

Lumped mass matrix:
```
M_lumped = ρAL/2 × diag([1, 1, L²/12, 1, 1, L²/12])
```

#### E.1.5 Implementation Code Template

```python
def ke_euler_bernoulli_2d(xe, ge):
    """
    Compute 6x6 stiffness matrix for Euler-Bernoulli beam.
    
    Parameters
    ----------
    xe : ndarray, shape (2, 2)
        Node coordinates [[x1, y1], [x2, y2]]
    ge : ndarray
        Material properties [E, A, I]
        E = Young's modulus
        A = Cross-section area
        I = Moment of inertia
    
    Returns
    -------
    ke : ndarray, shape (6, 6)
        Element stiffness matrix in global coordinates
    """
    # Extract properties
    E, A, I = ge[0], ge[1], ge[2]
    
    # Compute length and direction cosines
    dx = xe[1, 0] - xe[0, 0]
    dy = xe[1, 1] - xe[0, 1]
    L = np.sqrt(dx**2 + dy**2)
    c = dx / L  # cos(θ)
    s = dy / L  # sin(θ)
    
    # Local stiffness matrix
    k_local = np.zeros((6, 6))
    
    # Axial terms
    k_local[0, 0] = E * A / L
    k_local[0, 3] = -E * A / L
    k_local[3, 0] = -E * A / L
    k_local[3, 3] = E * A / L
    
    # Bending terms
    EI_L3 = E * I / L**3
    EI_L2 = E * I / L**2
    EI_L = E * I / L
    
    k_local[1, 1] = 12 * EI_L3
    k_local[1, 2] = 6 * EI_L2
    k_local[1, 4] = -12 * EI_L3
    k_local[1, 5] = 6 * EI_L2
    
    k_local[2, 1] = 6 * EI_L2
    k_local[2, 2] = 4 * EI_L
    k_local[2, 4] = -6 * EI_L2
    k_local[2, 5] = 2 * EI_L
    
    k_local[4, 1] = -12 * EI_L3
    k_local[4, 2] = -6 * EI_L2
    k_local[4, 4] = 12 * EI_L3
    k_local[4, 5] = -6 * EI_L2
    
    k_local[5, 1] = 6 * EI_L2
    k_local[5, 2] = 2 * EI_L
    k_local[5, 4] = -6 * EI_L2
    k_local[5, 5] = 4 * EI_L
    
    # Transformation matrix
    T = np.array([
        [c,  s, 0,  0,  0, 0],
        [-s, c, 0,  0,  0, 0],
        [0,  0, 1,  0,  0, 0],
        [0,  0, 0,  c,  s, 0],
        [0,  0, 0, -s,  c, 0],
        [0,  0, 0,  0,  0, 1]
    ])
    
    # Transform to global
    ke = T.T @ k_local @ T
    
    return ke
```

### E.2 Timoshenko Beam Element (2D)

#### E.2.1 Kinematics

The Timoshenko beam includes shear deformation:
```
γ = dv/dx - θ
```

Where:
- γ = shear strain
- θ = rotation (independent of slope)

#### E.2.2 Shear Correction Factor

For rectangular cross-section:
```
κ = 5/6 ≈ 0.833
```

For circular cross-section:
```
κ = 6(1 + ν) / (7 + 6ν)
```

For I-section:
```
κ ≈ A_web / A_total
```

#### E.2.3 Stiffness Matrix (with Shear Deformation)

Define:
```
Φ = 12EI / (κGAL²)
```

Stiffness coefficients:
```
K_bb = 12EI / [L³(1 + Φ)]
K_bc = 6EI / [L²(1 + Φ)]
K_cc = (4 + Φ)EI / [L(1 + Φ)]
K_ce = (2 - Φ)EI / [L(1 + Φ)]
```

#### E.2.4 Reduced Integration

For thin beams, Φ → 0 and Timoshenko → Euler-Bernoulli.
For thick beams (L/h < 10), Timoshenko gives more accurate results.

Use 1-point reduced integration to avoid shear locking:
```
∫₋₁¹ f(ξ) dξ ≈ 2 × f(0)
```

### E.3 Shell Element (MITC4)

#### E.3.1 Geometry

The MITC4 shell element uses a bilinear isoparametric formulation with:
- 4 nodes
- 5 DOFs per node: (u, v, w, θₓ, θᵧ) or (u, v, w, α, β)
- Total 20 DOFs (or 24 with drilling DOF)

#### E.3.2 Coordinate Systems

Three coordinate systems:
1. **Global**: (X, Y, Z) - fixed reference
2. **Local**: (x, y, z) - element-specific, z normal to mid-surface
3. **Natural**: (ξ, η) - -1 ≤ ξ, η ≤ 1

#### E.3.3 Shape Functions

Bilinear shape functions:
```
N₁(ξ, η) = (1 - ξ)(1 - η) / 4
N₂(ξ, η) = (1 + ξ)(1 - η) / 4
N₃(ξ, η) = (1 + ξ)(1 + η) / 4
N₄(ξ, η) = (1 - ξ)(1 + η) / 4
```

#### E.3.4 MITC Assumed Strain Interpolation

To avoid transverse shear locking, MITC4 uses assumed strains:

Sampling points:
```
Point A: (ξ = 0, η = -1)  - bottom edge midpoint
Point B: (ξ = 1, η = 0)   - right edge midpoint
Point C: (ξ = 0, η = 1)   - top edge midpoint
Point D: (ξ = -1, η = 0)  - left edge midpoint
```

Assumed shear strains:
```
γₓₖ = (1 - η)/2 × γₓₖ|_A + (1 + η)/2 × γₓₖ|_C
γᵧₖ = (1 - ξ)/2 × γᵧₖ|_D + (1 + ξ)/2 × γᵧₖ|_B
```

#### E.3.5 Stiffness Matrix Decomposition

Shell stiffness combines:
```
K = K_membrane + K_bending + K_shear
```

Membrane stiffness (constant stress triangle assembly):
```
K_membrane = ∫∫ Bₘᵀ Dₘ Bₘ t dA
```

Bending stiffness:
```
K_bending = ∫∫ Bᵦᵀ Dᵦ Bᵦ t³/12 dA
```

Shear stiffness (with MITC):
```
K_shear = ∫∫ B̄ₛᵀ Dₛ B̄ₛ κt dA
```

Where B̄ₛ uses the assumed shear strain interpolation.

### E.4 DKT Shell Triangle

#### E.4.1 Overview

Discrete Kirchhoff Triangle (DKT):
- 3 nodes, 3 DOFs per node (w, θₓ, θᵧ)
- 9 DOFs total
- Kirchhoff hypothesis enforced discretely at selected points

#### E.4.2 Displacement Field

Transverse displacement:
```
w(ξ, η) = Σᵢ Nᵢ wᵢ + Σₖ Nₖ wₖ
```

Where:
- Nᵢ = corner node shape functions
- Nₖ = edge midpoint shape functions (eliminated)

#### E.4.3 Rotation Interpolation

Rotations interpolated using:
```
βₓ(ξ, η) = Σⱼ Hⱼₓ(ξ, η) × qⱼ
βᵧ(ξ, η) = Σⱼ Hⱼᵧ(ξ, η) × qⱼ
```

Where qⱼ = [w₁, θₓ₁, θᵧ₁, w₂, θₓ₂, θᵧ₂, w₃, θₓ₃, θᵧ₃]ᵀ

#### E.4.4 Kirchhoff Constraint

At edge midpoints (k = 4, 5, 6):
```
∂w/∂s|ₖ = -βₛ|ₖ
```

Where s = edge tangent direction.

This constraint relates edge midpoint DOFs to corner DOFs, reducing 15 DOFs to 9.

### E.5 Higher-Order Elements

#### E.5.1 Six-Node Triangle (T6)

Shape functions (quadratic):
```
N₁ = (2L₁ - 1)L₁
N₂ = (2L₂ - 1)L₂
N₃ = (2L₃ - 1)L₃
N₄ = 4L₁L₂
N₅ = 4L₂L₃
N₆ = 4L₃L₁
```

Where L₁, L₂, L₃ are area coordinates with L₁ + L₂ + L₃ = 1.

Integration: 3-point rule for full integration, 1-point for reduced.

#### E.5.2 Eight-Node Quadrilateral (Q8)

Serendipity shape functions:
```
Corner nodes (i = 1, 2, 3, 4):
Nᵢ = (1 + ξξᵢ)(1 + ηηᵢ)(ξξᵢ + ηηᵢ - 1) / 4

Midside nodes (i = 5, 7):
Nᵢ = (1 - ξ²)(1 + ηηᵢ) / 2

Midside nodes (i = 6, 8):
Nᵢ = (1 + ξξᵢ)(1 - η²) / 2
```

Integration: 3×3 for full, 2×2 for reduced.

#### E.5.3 Ten-Node Tetrahedron (T10)

Volume coordinates: L₁, L₂, L₃, L₄ with L₁ + L₂ + L₃ + L₄ = 1.

Shape functions:
```
Corner nodes (i = 1, 2, 3, 4):
Nᵢ = (2Lᵢ - 1)Lᵢ

Edge midpoint nodes (i = 5, ..., 10):
Nᵢⱼ = 4LᵢLⱼ  (for edge connecting corners i and j)
```

Integration: 4-point rule for full integration.

#### E.5.4 Twenty-Node Hexahedron (H20)

Serendipity shape functions in 3D:
```
Corner nodes (ξᵢ, ηᵢ, ζᵢ = ±1):
Nᵢ = (1 + ξξᵢ)(1 + ηηᵢ)(1 + ζζᵢ)(ξξᵢ + ηηᵢ + ζζᵢ - 2) / 8

Edge midpoint nodes:
Various forms depending on which coordinate is zero
```

Integration: 3×3×3 for full, 2×2×2 for reduced.

---

## Appendix F: Material Model Mathematical Formulations

### F.1 Steel02 (Giuffré-Menegotto-Pinto Model)

#### F.1.1 Backbone Curve

The stress-strain relationship:
```
σ = b × ε + (1 - b) × ε / [1 + |ε|ᴿ]^(1/R)
```

Where:
- b = E_hard / E₀ (ratio of hardening to initial stiffness)
- R = curvature parameter (controls transition sharpness)
- ε = (ε - ε_r) / (ε₀ - ε_r) (normalized strain)
- σ = (σ - σ_r) / (σ₀ - σ_r) (normalized stress)

#### F.1.2 Bauschinger Effect

The R parameter evolves with plastic strain:
```
R = R₀ - (a₁ × ξ) / (a₂ + ξ)
```

Where:
- R₀ = initial R value (typ. 15-20)
- ξ = cumulative plastic strain
- a₁, a₂ = calibration parameters

#### F.1.3 Isotropic Hardening

Optional isotropic hardening shifts the asymptotes:
```
σ_shift = a₃ × [σ_y × (a₄ × ε_p_max + ε_p)]
```

#### F.1.4 Implementation Algorithm

```python
def stress_steel02(eps, state, props):
    """
    Steel02 material state determination.
    
    Parameters
    ----------
    eps : float
        Current strain
    state : dict
        Material state variables:
        - eps_r, sig_r: reversal point
        - eps_0, sig_0: target point
        - eps_min, eps_max: strain history
        - kon: loading indicator
    props : dict
        Material properties:
        - Fy: yield stress
        - E0: initial modulus
        - b: hardening ratio
        - R0, cR1, cR2: Bauschinger parameters
        - a1, a2, a3, a4: isotropic hardening
    
    Returns
    -------
    sig : float
        Current stress
    Et : float
        Tangent modulus
    state : dict
        Updated state
    """
    # Unpack properties
    Fy = props['Fy']
    E0 = props['E0']
    b = props['b']
    R0, cR1, cR2 = props['R0'], props['cR1'], props['cR2']
    
    # Unpack state
    eps_r = state['eps_r']
    sig_r = state['sig_r']
    eps_0 = state['eps_0']
    sig_0 = state['sig_0']
    
    # Compute shifted yield
    epsyP = Fy / E0
    
    # Check for reversal
    deps = eps - state['eps_prev']
    if state['kon'] == 0:  # initial loading
        # ... initialize loading direction
        pass
    elif state['kon'] * deps < 0:  # reversal
        # Update reversal point
        eps_r = state['eps_prev']
        sig_r = state['sig_prev']
        # Compute new target
        # ... 
        
    # Compute R with Bauschinger
    xi = abs((state['eps_max'] - state['eps_min']) / (2 * epsyP))
    R = R0 - cR1 * xi / (cR2 + xi)
    
    # Normalized strain
    eps_norm = (eps - eps_r) / (eps_0 - eps_r)
    
    # Menegotto-Pinto curve
    sig_norm = b * eps_norm + (1 - b) * eps_norm / (1 + abs(eps_norm)**R)**(1/R)
    
    # Denormalize
    sig = sig_r + sig_norm * (sig_0 - sig_r)
    
    # Tangent modulus
    Et = b * E0 + (1 - b) * E0 / (1 + abs(eps_norm)**R)**(1 + 1/R)
    
    # Update state
    state['eps_prev'] = eps
    state['sig_prev'] = sig
    state['eps_max'] = max(state['eps_max'], eps)
    state['eps_min'] = min(state['eps_min'], eps)
    
    return sig, Et, state
```

### F.2 Concrete02 Model

#### F.2.1 Compression Envelope

Modified Kent-Park model:
```
For ε < ε_c0:
    σ_c = f_c × [2(ε/ε_c0) - (ε/ε_c0)²]

For ε_c0 ≤ ε < ε_cu:
    σ_c = f_c × [1 - Z(ε - ε_c0)]

For ε ≥ ε_cu:
    σ_c = 0.2 × f_c  (residual)
```

Where:
- f_c = compressive strength (negative)
- ε_c0 = strain at peak (typ. -0.002)
- Z = softening slope = 0.5 / (ε_50u - ε_c0)
- ε_50u = strain at 50% strength

#### F.2.2 Tension Envelope

Linear softening:
```
For 0 < ε < ε_t0:
    σ_t = E_t × ε

For ε_t0 ≤ ε < ε_tu:
    σ_t = f_t × [1 - (ε - ε_t0) / (ε_tu - ε_t0)]

For ε ≥ ε_tu:
    σ_t = 0
```

Where:
- f_t = tensile strength (typ. 0.1 × f_c)
- ε_t0 = cracking strain = f_t / E₀
- ε_tu = ultimate tensile strain

#### F.2.3 Unloading/Reloading Rules

Compression unloading:
```
σ = σ_un + E_un × (ε - ε_un)
E_un = f_c / ε_c0  (initial unloading)
E_un = σ_un / (ε_un - ε_pl)  (subsequent)
```

Where:
- ε_pl = plastic strain
- σ_un, ε_un = unloading point

Tension-compression transition:
```
Reloading to compression after tension crack closes at ε = 0
```

### F.3 Fiber Section Integration

#### F.3.1 Section Force-Deformation

For a 2D beam section with deformation vector φ = [ε₀, κ]:
```
ε(y) = ε₀ - y × κ
```

Section forces:
```
N = ∫_A σ dA = Σᵢ Aᵢ × σᵢ
M = ∫_A -y × σ dA = Σᵢ -yᵢ × Aᵢ × σᵢ
```

#### F.3.2 Section Tangent Stiffness

```
k_section = [∂N/∂ε₀   ∂N/∂κ ]   [Σ Aᵢ Eᵢ     -Σ yᵢ Aᵢ Eᵢ]
            [∂M/∂ε₀   ∂M/∂κ ] = [-Σ yᵢ Aᵢ Eᵢ  Σ yᵢ² Aᵢ Eᵢ]
```

Where Eᵢ = tangent modulus of fiber i.

#### F.3.3 Implementation

```python
def fiber_section_response(phi, fibers):
    """
    Compute fiber section response.
    
    Parameters
    ----------
    phi : ndarray, shape (2,)
        Section deformations [eps_0, kappa]
    fibers : list of dict
        Each fiber: {'y': y_coord, 'A': area, 'material': mat_obj}
    
    Returns
    -------
    s : ndarray, shape (2,)
        Section forces [N, M]
    ks : ndarray, shape (2, 2)
        Section tangent stiffness
    """
    eps_0, kappa = phi
    
    N = 0.0
    M = 0.0
    k11 = 0.0  # dN/d_eps0
    k12 = 0.0  # dN/d_kappa
    k22 = 0.0  # dM/d_kappa
    
    for fiber in fibers:
        y = fiber['y']
        A = fiber['A']
        mat = fiber['material']
        
        # Fiber strain
        eps = eps_0 - y * kappa
        
        # Material response
        sigma, Et = mat.get_stress(eps)
        
        # Section forces
        N += A * sigma
        M += -y * A * sigma
        
        # Tangent contributions
        k11 += A * Et
        k12 += -y * A * Et
        k22 += y**2 * A * Et
    
    s = np.array([N, M])
    ks = np.array([[k11, k12],
                   [k12, k22]])
    
    return s, ks
```

---

## Appendix G: Time Integration Algorithm Details

### G.1 Newmark-β Method (Implicit)

#### G.1.1 Algorithm Steps

Given: M, C, K, p(t), u_n, v_n, a_n, dt

**Step 1: Form effective stiffness**
```
K̂ = K + (γ / β / dt) × C + (1 / β / dt²) × M
```

**Step 2: Form effective load**
```
p̂ = p_{n+1} + M × [(1/β/dt²) × u_n + (1/β/dt) × v_n + (1/2β - 1) × a_n]
           + C × [(γ/β/dt) × u_n + (γ/β - 1) × v_n + dt × (γ/2β - 1) × a_n]
```

**Step 3: Solve for displacement**
```
K̂ × u_{n+1} = p̂
```

**Step 4: Update velocity and acceleration**
```
a_{n+1} = (1/β/dt²) × (u_{n+1} - u_n) - (1/β/dt) × v_n - (1/2β - 1) × a_n
v_{n+1} = v_n + dt × [(1 - γ) × a_n + γ × a_{n+1}]
```

#### G.1.2 Parameters

| Method | β | γ | Order | Stability |
|--------|---|---|-------|-----------|
| Average acceleration | 1/4 | 1/2 | 2 | Unconditional |
| Linear acceleration | 1/6 | 1/2 | 2 | Conditional |
| Fox-Goodwin | 1/12 | 1/2 | 4* | Conditional |
| Central difference | 0 | 1/2 | 2 | Conditional |

*4th order for undamped SDOF

#### G.1.3 Stability Condition

For conditional methods:
```
dt ≤ dt_crit = T_min / π × √(γ/2 - β) / (γ/2)
```

For average acceleration (unconditional):
```
dt ≤ T_min / 10  (accuracy recommendation)
```

### G.2 HHT-α Method

#### G.2.1 Modified Equations

The HHT method evaluates equilibrium at time (1 + α)dt:
```
M × a_{n+1} + (1 + α) × C × v_{n+1} - α × C × v_n 
+ (1 + α) × K × u_{n+1} - α × K × u_n = (1 + α) × p_{n+1} - α × p_n
```

#### G.2.2 Algorithm

Same as Newmark with modified effective stiffness:
```
K̂ = (1 + α) × K + (γ / β / dt) × (1 + α) × C + (1 / β / dt²) × M
```

Modified effective load includes α × (previous step forces).

#### G.2.3 Parameters

```
-1/3 ≤ α ≤ 0  (α = 0 recovers Newmark)
β = (1 - α)² / 4
γ = (1 - 2α) / 2
```

Typical choice: α = -0.05 (slight numerical damping)

### G.3 Central Difference (Explicit)

#### G.3.1 Algorithm

**Step 1: Compute acceleration**
```
a_n = M⁻¹ × (p_n - C × v_n - K × u_n)
```

**Step 2: Update velocity (half step)**
```
v_{n+1/2} = v_{n-1/2} + dt × a_n
```

**Step 3: Update displacement**
```
u_{n+1} = u_n + dt × v_{n+1/2}
```

#### G.3.2 Lumped Mass

For explicit methods, lumped (diagonal) mass is essential:
```
M_ii = Σ_e M_e × row_sum_contribution
```

No matrix factorization needed, only vector operations.

#### G.3.3 Stability

CFL condition:
```
dt ≤ dt_crit = 2 / ω_max = L_min / c
```

Where:
- ω_max = maximum natural frequency
- L_min = smallest element dimension
- c = wave speed = √(E/ρ)

### G.4 Arc-Length Method

#### G.4.1 Spherical Constraint

Constrain solution to lie on hypersphere:
```
Δu² + ψ² × Δλ² × |p_ref|² = Δl²
```

Where:
- Δu = displacement increment
- Δλ = load factor increment
- ψ = scaling parameter
- Δl = arc length

#### G.4.2 Predictor Step

Initial guess:
```
Δu_pred = ±Δl × û / √(û² + ψ² × |p_ref|²)
Δλ_pred = ±Δl × 1 / √(û² / |p_ref|² + ψ²)
```

Sign determined by dot product with previous increment.

#### G.4.3 Corrector Iterations

Newton-Raphson with constraint:
```
[K  -p_ref] [δu ]   [r - f(u)]
[aᵀ  b    ] [δλ ] = [c       ]
```

Where:
- a = ∂g/∂u = 2Δu
- b = ∂g/∂λ = 2ψ²Δλ|p_ref|²
- c = Δl² - Δu² - ψ²Δλ²|p_ref|²

---

## Appendix H: Validation Test Cases

### H.1 NAFEMS Benchmarks

#### H.1.1 Linear Elastic Benchmarks

| Test ID | Description | Element Types | Reference |
|---------|-------------|---------------|-----------|
| LE1 | Elliptic membrane | Q4, Q8, T3, T6 | NAFEMS R0016 |
| LE2 | Cylindrical shell patch | Shell | NAFEMS R0016 |
| LE3 | Hemisphere-point loads | Shell | NAFEMS R0016 |
| LE5 | Z-section cantilever | Shell | NAFEMS R0016 |
| LE6 | Skew plate | Shell | NAFEMS R0016 |
| LE10 | Thick cylinder | 3D solid | NAFEMS R0016 |
| LE11 | Solid cylinder | 3D solid | NAFEMS R0016 |

#### H.1.2 Natural Frequency Benchmarks

| Test ID | Description | Modes | Reference |
|---------|-------------|-------|-----------|
| FV2 | Pin-ended cross beam | 1-12 | NAFEMS FFV |
| FV4 | Free-free cross beam | 1-12 | NAFEMS FFV |
| FV12 | Free square thin plate | 1-6 | NAFEMS FFV |
| FV22 | Clamped thick cylinder | 1-20 | NAFEMS FFV |
| FV32 | Cantilevered tapered membrane | 1-12 | NAFEMS FFV |
| FV41 | Free thin curved shell | 1-12 | NAFEMS FFV |
| FV52 | Solid cylinder | 1-20 | NAFEMS FFV |

#### H.1.3 Nonlinear Benchmarks

| Test ID | Description | Type | Reference |
|---------|-------------|------|-----------|
| NL1 | Z-section beam | Geometric NL | NAFEMS NL1 |
| NL2 | Clamped hemispherical dome | Geometric NL | NAFEMS NL2 |
| NL3 | Hemispherical shell | Large rotation | NAFEMS NL3 |
| NL4 | Snap-through of arc | Instability | NAFEMS NL4 |
| NL5 | Two-bar truss snap-through | Geometric NL | Custom |
| NL6 | Lee frame | Combined NL | Custom |

### H.2 OpenSees Verification Examples

#### H.2.1 Static Examples
```python
# Example: Cantilever beam
def test_cantilever_opensees():
    # Build model in OpenSeesPy
    import openseespy.opensees as ops
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    # ... create nodes, elements, BCs, loads
    ops.analyze(1)
    u_opensees = ops.nodeDisp(2, 2)
    
    # Build model in femlabpy
    # ... same geometry
    result = fp.elastic(T, X, G, C, P, dof=3)
    u_femlabpy = result['u'][1, 1]
    
    assert abs(u_femlabpy - u_opensees) / abs(u_opensees) < 0.01
```

#### H.2.2 Dynamic Examples
```python
# Example: SDOF oscillator
def test_sdof_newmark_opensees():
    # Natural period T = 1.0 s, damping 5%
    # Harmonic load p(t) = sin(2π/T * t)
    # Compare time history over 10 cycles
    
    # OpenSeesPy
    ops.timeSeries('Trig', 1, 0.0, 10.0, 1.0)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0, 0.0, 0.0)
    # ... analysis
    
    # femlabpy
    result = fp.solve_newmark(M, C, K, p_func, u0, v0, dt, nsteps)
    
    # Compare peak displacement
    assert abs(result.u.max() - u_opensees_max) < 0.01
```

### H.3 CalculiX Verification Examples

#### H.3.1 Static Examples
```
*NODE
1, 0.0, 0.0
2, 1.0, 0.0
...
*ELEMENT, TYPE=CPS4
1, 1, 2, 6, 5
...
*MATERIAL, NAME=STEEL
*ELASTIC
210000.0, 0.3
*STEP
*STATIC
*BOUNDARY
1, 1, 2
...
*CLOAD
10, 1, 1000.0
*END STEP
```

Compare displacement results from CCX .frd file with femlabpy.

---

## Appendix I: Code Quality Standards

### I.1 Python Style Guide

Follow PEP 8 with these project-specific conventions:

1. **Docstrings**: NumPy style with Mathematical Formulation and Algorithm sections
2. **Type hints**: Use for all public functions
3. **Line length**: 88 characters (Black formatter)
4. **Imports**: Group by standard library, third-party, local

### I.2 Testing Standards

1. **Coverage**: Minimum 80% line coverage
2. **Unit tests**: One test file per module
3. **Integration tests**: Cross-module workflows
4. **Benchmark tests**: Against reference solutions

### I.3 Documentation Standards

1. **API reference**: Auto-generated from docstrings
2. **User guide**: Tutorial-style with examples
3. **Theory manual**: Mathematical derivations
4. **Developer guide**: Architecture and contribution

### I.4 Performance Standards

1. **Assembly**: Use vectorized operations (np.einsum)
2. **Sparse matrices**: Use scipy.sparse for large problems
3. **Solve**: Auto-select direct vs iterative based on size
4. **Memory**: Avoid unnecessary copies

---

## Appendix J: Dependency Matrix

### J.1 Task Dependencies (Graph)

```
Phase 1.1 (Beam 2D) ──┐
                      ├──> Phase 1.3 (Shell)
Phase 1.2 (Beam 3D) ──┘

Phase 1.4 (Spring) ──> Phase 2.3 (Hysteretic)

Phase 2.1 (Steel) ──┐
                    ├──> Phase 2.4 (Fiber Section)
Phase 2.2 (Concrete) ─┘

Phase 3.1 (Geom NL) ──┐
                      ├──> Phase 3.2 (Arc-Length)
Phase 3.3 (Buckling) ──┘

Phase 3.4 (Spectrum) <── Phase 5.3 (Modal Superposition)
```

### J.2 Module Dependencies

```python
# femlabpy dependency graph
femlabpy/
├── _helpers.py         # No internal deps
├── types.py            # No internal deps
├── elements/
│   ├── triangles.py    # _helpers
│   ├── quads.py        # _helpers
│   ├── bars.py         # _helpers
│   ├── solids.py       # _helpers
│   ├── beams2d.py      # _helpers (NEW)
│   ├── beams3d.py      # _helpers (NEW)
│   └── shells.py       # _helpers (NEW)
├── materials/
│   ├── invariants.py   # _helpers
│   ├── plasticity.py   # invariants
│   ├── steel.py        # _helpers (NEW)
│   ├── concrete.py     # _helpers (NEW)
│   ├── hysteretic.py   # _helpers (NEW)
│   └── sections.py     # steel, concrete (NEW)
├── assembly.py         # elements
├── boundary.py         # _helpers
├── loads.py            # _helpers
├── solvers.py          # elements, boundary, loads
├── modal.py            # _helpers
├── dynamics.py         # _helpers, modal
├── damping.py          # _helpers
├── periodic.py         # _helpers
└── analysis/
    ├── corotational.py # elements (NEW)
    ├── arclength.py    # solvers (NEW)
    ├── buckling.py     # assembly, modal (NEW)
    ├── spectrum.py     # modal (NEW)
    └── thermal.py      # elements (NEW)
```

---

## Appendix K: Complete OpenSees Source File Catalog

This appendix provides a complete listing of all source files in the OpenSees repository
that are relevant to femlabpy implementation. Each file is categorized and annotated
with implementation priority and complexity.

### K.1 Element Source Files

#### K.1.1 Beam-Column Elements

##### 2D Beam Elements
```
SRC/element/beam2d/
├── beam2d.h                          # Abstract base class for 2D beams
├── beam2d01.cpp                      # Basic 2D beam element (Priority: CRITICAL)
├── beam2d01.h                        # 
├── beam2d02.cpp                      # 2D beam with hinges (Priority: HIGH)
├── beam2d02.h                        # 
├── beam2d03.cpp                      # 2D beam with end releases (Priority: HIGH)
├── beam2d03.h                        # 
├── beam2d04.cpp                      # 2D beam with nonlinear springs (Priority: MEDIUM)
└── beam2d04.h                        # 

Key classes to implement:
- Beam2d01: Basic Euler-Bernoulli beam element
  - Methods: setDomain(), commitState(), revertToLastCommit()
  - Stiffness: getInitialStiff(), getTangentStiff()
  - Resisting force: getResistingForce(), getResistingForceIncInertia()
  - DOFs: 3 per node (u, v, θ)
  
- Beam2d02: Beam with rotational springs at ends
  - Additional parameters: spring stiffness k1, k2
  - Condensation of internal DOFs
```

##### 3D Beam Elements
```
SRC/element/beam3d/
├── beam3d.h                          # Abstract base class for 3D beams
├── beam3d01.cpp                      # Basic 3D beam (Priority: HIGH)
├── beam3d01.h                        #
├── beam3d02.cpp                      # 3D beam with hinges (Priority: MEDIUM)
└── beam3d02.h                        #

Key implementation notes:
- 3D beams require orientation definition (local axes)
- 6 DOFs per node: (ux, uy, uz, θx, θy, θz)
- Transformation matrix T is 12×12
- Geometric stiffness for buckling (P-δ, P-Δ effects)
```

##### Elastic Beam-Column Elements
```
SRC/element/elasticBeamColumn/
├── ElasticBeam2d.cpp                 # 2D elastic beam (Priority: CRITICAL)
├── ElasticBeam2d.h                   #
├── ElasticBeam3d.cpp                 # 3D elastic beam (Priority: HIGH)
├── ElasticBeam3d.h                   #
├── ElasticTimoshenkoBeam2d.cpp       # 2D Timoshenko (Priority: HIGH)
├── ElasticTimoshenkoBeam2d.h         #
├── ElasticTimoshenkoBeam3d.cpp       # 3D Timoshenko (Priority: MEDIUM)
├── ElasticTimoshenkoBeam3d.h         #
├── ModElasticBeam2d.cpp              # Modified elastic beam (Priority: LOW)
├── ModElasticBeam2d.h                #
├── ModElasticBeam3d.cpp              # (Priority: LOW)
└── ModElasticBeam3d.h                #

Timoshenko beam implementation:
- Includes shear deformation: γ = G*A_s where A_s = κ*A
- Shear correction factor κ depends on cross-section shape
- Reduced integration to avoid shear locking
- Mass matrix includes rotary inertia
```

##### Displacement-Based Beam-Column Elements
```
SRC/element/dispBeamColumn/
├── DispBeamColumn2d.cpp              # Displacement-based formulation (Priority: MEDIUM)
├── DispBeamColumn2d.h                #
├── DispBeamColumn3d.cpp              # (Priority: MEDIUM)
├── DispBeamColumn3d.h                #
├── DispBeamColumn2dThermal.cpp       # With thermal effects (Priority: LOW)
├── DispBeamColumn2dThermal.h         #
├── DispBeamColumnNL2d.cpp            # Geometrically nonlinear (Priority: MEDIUM)
├── DispBeamColumnNL2d.h              #
├── TimoshenkoBeamColumn2d.cpp        # Timoshenko formulation (Priority: MEDIUM)
└── TimoshenkoBeamColumn2d.h          #

Displacement-based formulation:
- Assumes displacement field along element
- Uses Gauss-Lobatto integration points
- Distributed plasticity through fiber sections
- Less accurate than force-based for softening
```

##### Force-Based Beam-Column Elements
```
SRC/element/forceBeamColumn/
├── ForceBeamColumn2d.cpp             # Force-based formulation (Priority: HIGH)
├── ForceBeamColumn2d.h               #
├── ForceBeamColumn3d.cpp             # (Priority: HIGH)
├── ForceBeamColumn3d.h               #
├── ForceBeamColumnCBDI2d.cpp         # With curvature-based displacement (Priority: LOW)
├── ForceBeamColumnCBDI2d.h           #
├── ForceBeamColumnCBDI3d.cpp         # (Priority: LOW)
├── ForceBeamColumnCBDI3d.h           #
├── ForceBeamColumnWarping2d.cpp      # With warping (Priority: LOW)
├── ForceBeamColumnWarping2d.h        #
├── GradientInelasticBeamColumn2d.cpp # Gradient plasticity (Priority: LOW)
├── GradientInelasticBeamColumn2d.h   #
├── GradientInelasticBeamColumn3d.cpp # (Priority: LOW)
├── GradientInelasticBeamColumn3d.h   #
└── TclForceBeamColumnCommand.cpp     # Tcl interface

Force-based formulation:
- Interpolates force along element (exact equilibrium)
- Uses Gauss-Lobatto integration
- More accurate for strain softening
- Requires state determination iterative procedure
```

#### K.1.2 Shell Elements
```
SRC/element/shell/
├── ASDShellQ4.cpp                    # 4-node MITC shell (Priority: CRITICAL)
├── ASDShellQ4.h                      #
├── ShellDKGQ.cpp                     # DKG quadrilateral (Priority: HIGH)
├── ShellDKGQ.h                       #
├── ShellDKGT.cpp                     # DKG triangle (Priority: HIGH)
├── ShellDKGT.h                       #
├── ShellMITC4.cpp                    # MITC4 element (Priority: HIGH)
├── ShellMITC4.h                      #
├── ShellMITC4Thermal.cpp             # With thermal (Priority: LOW)
├── ShellMITC4Thermal.h               #
├── ShellMITC9.cpp                    # 9-node MITC (Priority: MEDIUM)
├── ShellMITC9.h                      #
├── ShellNLDKGQ.cpp                   # Nonlinear DKG quad (Priority: MEDIUM)
├── ShellNLDKGQ.h                     #
├── ShellNLDKGT.cpp                   # Nonlinear DKG tri (Priority: MEDIUM)
├── ShellNLDKGT.h                     #
├── R3vectors.cpp                     # Vector utilities
├── R3vectors.h                       #
├── ShellCorotational.cpp             # Corotational formulation (Priority: MEDIUM)
└── ShellCorotational.h               #

MITC4 implementation notes:
- Mixed Interpolation of Tensorial Components
- Eliminates transverse shear locking
- 4 nodes, 5 or 6 DOFs per node
- Assumed strain field for shear
- Tying points at edge midpoints

DKT/DKG implementation notes:
- Discrete Kirchhoff Triangle/Quadrilateral
- No transverse shear (thin shell assumption)
- 3 DOFs per node: (w, θx, θy)
- Edge rotations constrained to match deflection slope
```

#### K.1.3 Solid Elements
```
SRC/element/brick/
├── Brick.cpp                         # 8-node brick (Priority: EXISTS - H8)
├── Brick.h                           #
├── BbarBrick.cpp                     # B-bar formulation (Priority: HIGH)
├── BbarBrick.h                       #
├── Brick20N.cpp                      # 20-node brick (Priority: HIGH)
├── Brick20N.h                        #
├── Brick27N.cpp                      # 27-node brick (Priority: MEDIUM)
├── Brick27N.h                        #
├── BrickSurfaceLoad.cpp              # Surface load element (Priority: LOW)
├── BrickSurfaceLoad.h                #
├── BrickUP.cpp                       # Coupled U-P (Priority: LOW)
├── BrickUP.h                         #
├── TclBrickCommand.cpp               # Tcl interface
└── TwentyEightNodeBrickUP.cpp        # 28-node U-P brick

B-bar formulation:
- Selective reduced integration
- Avoids volumetric locking for incompressible materials
- Deviatoric strain: ε_dev = ε - (1/3)tr(ε)I
- Volumetric strain computed at center: ε_vol = (1/3)tr(ε)|_center

SRC/element/tetrahedron/
├── FourNodeTetrahedron.cpp           # 4-node tet (Priority: EXISTS - T4)
├── FourNodeTetrahedron.h             #
├── TenNodeTetrahedron.cpp            # 10-node tet (Priority: HIGH)
└── TenNodeTetrahedron.h              #

10-node tetrahedron:
- Quadratic displacement field
- 4 integration points (4-point rule)
- Better accuracy than T4 for bending
```

#### K.1.4 Plane Elements
```
SRC/element/fourNodeQuad/
├── FourNodeQuad.cpp                  # 4-node quad (Priority: EXISTS - Q4)
├── FourNodeQuad.h                    #
├── EnhancedQuad.cpp                  # Enhanced strain quad (Priority: HIGH)
├── EnhancedQuad.h                    #
├── FourNodeQuadUP.cpp                # Coupled U-P (Priority: LOW)
├── FourNodeQuadUP.h                  #
├── NineNodeMixedQuad.cpp             # 9-node mixed (Priority: MEDIUM)
├── NineNodeMixedQuad.h               #
├── NineNodeQuad.cpp                  # 9-node quad (Priority: MEDIUM)
├── NineNodeQuad.h                    #
├── NineNodeQuadUP.cpp                # 9-node U-P (Priority: LOW)
└── NineNodeQuadUP.h                  #

Enhanced strain formulation:
- Adds incompatible modes to displacement field
- Eliminates shear and volumetric locking
- Internal DOFs (α modes) condensed out
- Passes patch test

SRC/element/tri31/
├── Tri31.cpp                         # 3-node triangle (Priority: EXISTS - T3)
└── Tri31.h                           #

SRC/element/triangle/
├── Tri6N.cpp                         # 6-node triangle (Priority: HIGH)
└── Tri6N.h                           #
```

#### K.1.5 Truss and Bar Elements
```
SRC/element/truss/
├── CorotTruss.cpp                    # Corotational truss (Priority: HIGH)
├── CorotTruss.h                      #
├── CorotTrussSection.cpp             # With section (Priority: MEDIUM)
├── CorotTrussSection.h               #
├── Truss.cpp                         # Basic truss (Priority: EXISTS - Bar)
├── Truss.h                           #
├── Truss2.cpp                        # Alternative formulation (Priority: LOW)
├── Truss2.h                          #
├── TrussSection.cpp                  # With fiber section (Priority: MEDIUM)
└── TrussSection.h                    #

Corotational truss:
- Large displacement/rotation capability
- Rigid body rotation extracted
- Small strain assumption maintained
- Geometric stiffness included
```

#### K.1.6 Zero-Length Elements
```
SRC/element/zeroLength/
├── ZeroLength.cpp                    # General zero-length (Priority: HIGH)
├── ZeroLength.h                      #
├── ZeroLengthContact2D.cpp           # 2D contact (Priority: MEDIUM)
├── ZeroLengthContact2D.h             #
├── ZeroLengthContact3D.cpp           # 3D contact (Priority: MEDIUM)
├── ZeroLengthContact3D.h             #
├── ZeroLengthContactASDimplex.cpp    # Simplex contact (Priority: LOW)
├── ZeroLengthContactASDimplex.h      #
├── ZeroLengthContactNTS2D.cpp        # Node-to-segment (Priority: MEDIUM)
├── ZeroLengthContactNTS2D.h          #
├── ZeroLengthInterface2D.cpp         # Interface element (Priority: MEDIUM)
├── ZeroLengthInterface2D.h           #
├── ZeroLengthND.cpp                  # nD material (Priority: HIGH)
├── ZeroLengthND.h                    #
├── ZeroLengthRocking.cpp             # Rocking element (Priority: MEDIUM)
├── ZeroLengthRocking.h               #
├── ZeroLengthSection.cpp             # With section (Priority: HIGH)
└── ZeroLengthSection.h               #

Zero-length element uses:
- Springs, dampers, friction
- Concentrated plasticity hinges
- Gap/contact elements
- Structural connections
```

#### K.1.7 Joint Elements
```
SRC/element/joint/
├── BeamColumnJoint2d.cpp             # 2D beam-column joint (Priority: HIGH)
├── BeamColumnJoint2d.h               #
├── BeamColumnJoint3d.cpp             # 3D joint (Priority: MEDIUM)
├── BeamColumnJoint3d.h               #
├── Joint2D.cpp                       # Simplified 2D joint (Priority: HIGH)
├── Joint2D.h                         #
├── Joint3D.cpp                       # Simplified 3D joint (Priority: MEDIUM)
├── Joint3D.h                         #
├── MP_Joint2D.cpp                    # Multi-point joint (Priority: LOW)
└── MP_Joint2D.h                      #

Joint element features:
- Panel zone deformation
- Shear strength degradation
- Bond-slip modeling
- Rigid offset capability
```

#### K.1.8 Bearing Elements
```
SRC/element/elastomericBearing/
├── ElastomericBearing2d.cpp          # 2D elastomeric bearing (Priority: HIGH)
├── ElastomericBearing2d.h            #
├── ElastomericBearing3d.cpp          # 3D bearing (Priority: MEDIUM)
├── ElastomericBearing3d.h            #
├── ElastomericBearingBoucWen2d.cpp   # Bouc-Wen model (Priority: HIGH)
├── ElastomericBearingBoucWen2d.h     #
├── ElastomericBearingBoucWen3d.cpp   # (Priority: MEDIUM)
├── ElastomericBearingBoucWen3d.h     #
├── ElastomericBearingPlasticity2d.cpp # Plasticity model (Priority: HIGH)
├── ElastomericBearingPlasticity2d.h   #
├── ElastomericBearingPlasticity3d.cpp # (Priority: MEDIUM)
├── ElastomericBearingPlasticity3d.h   #
├── ElastomericBearingUFRP2d.cpp       # UFRP bearing (Priority: LOW)
└── ElastomericBearingUFRP2d.h         #

SRC/element/frictionBearing/
├── FlatSliderSimple2d.cpp             # Flat slider 2D (Priority: HIGH)
├── FlatSliderSimple2d.h               #
├── FlatSliderSimple3d.cpp             # Flat slider 3D (Priority: MEDIUM)
├── FlatSliderSimple3d.h               #
├── FrictionModel.cpp                  # Base friction model (Priority: HIGH)
├── FrictionModel.h                    #
├── Coulomb.cpp                        # Coulomb friction (Priority: HIGH)
├── Coulomb.h                          #
├── VelDependent.cpp                   # Velocity dependent (Priority: HIGH)
├── VelDependent.h                     #
├── VelDepMultiLinear.cpp              # Multi-linear (Priority: MEDIUM)
├── VelDepMultiLinear.h                #
├── VelNormalFrcDep.cpp                # Normal force dependent (Priority: MEDIUM)
├── VelNormalFrcDep.h                  #
├── SingleFPSimple2d.cpp               # Single FP 2D (Priority: HIGH)
├── SingleFPSimple2d.h                 #
├── SingleFPSimple3d.cpp               # Single FP 3D (Priority: MEDIUM)
├── SingleFPSimple3d.h                 #
├── MultiFP2d.cpp                      # Multi-surface FP (Priority: MEDIUM)
├── MultiFP2d.h                        #
├── TripleFrictionPendulum.cpp         # Triple FP (Priority: MEDIUM)
└── TripleFrictionPendulum.h           #

Bearing element features:
- Bidirectional interaction
- Vertical load variation
- Temperature effects
- Large displacement capability
```

### K.2 Material Source Files

#### K.2.1 Uniaxial Steel Materials
```
SRC/material/uniaxial/
├── Steel01.cpp                        # Bilinear steel (Priority: CRITICAL)
├── Steel01.h                          #
├── Steel02.cpp                        # Giuffré-Menegotto-Pinto (Priority: CRITICAL)
├── Steel02.h                          #
├── Steel02Fatigue.cpp                 # With fatigue (Priority: MEDIUM)
├── Steel02Fatigue.h                   #
├── Steel4.cpp                         # Comprehensive steel (Priority: HIGH)
├── Steel4.h                           #
├── SteelMPF.cpp                       # Modified MPF (Priority: MEDIUM)
├── SteelMPF.h                         #
├── SteelBRB.cpp                       # Buckling-restrained brace (Priority: HIGH)
├── SteelBRB.h                         #
├── DoddRestrepo.cpp                   # Dodd-Restrepo model (Priority: MEDIUM)
├── DoddRestrepo.h                     #
├── RambergOsgoodSteel.cpp             # Ramberg-Osgood (Priority: MEDIUM)
├── RambergOsgoodSteel.h               #
├── ReinforcingSteel.cpp               # Reinforcing bar (Priority: HIGH)
├── ReinforcingSteel.h                 #
├── UVCuniaxial.cpp                    # Updated Voce-Chaboche (Priority: MEDIUM)
└── UVCuniaxial.h                      #

Steel01 (bilinear with kinematic hardening):
- Parameters: Fy, E, b, a1, a2, a3, a4
- Kinematic hardening ratio b
- Optional isotropic hardening
- Asymmetric tension/compression

Steel02 (Menegotto-Pinto):
- Parameters: Fy, E, b, R0, cR1, cR2, a1, a2, a3, a4
- Smooth curved transition
- Bauschinger effect via R parameter
- Strain history tracking
```

#### K.2.2 Uniaxial Concrete Materials
```
SRC/material/uniaxial/
├── Concrete01.cpp                     # Kent-Park-Scott (Priority: CRITICAL)
├── Concrete01.h                       #
├── Concrete02.cpp                     # Linear tension softening (Priority: CRITICAL)
├── Concrete02.h                       #
├── Concrete04.cpp                     # Popovics (Priority: HIGH)
├── Concrete04.h                       #
├── Concrete06.cpp                     # Thorenfeldt (Priority: HIGH)
├── Concrete06.h                       #
├── Concrete07.cpp                     # Chang-Mander (Priority: MEDIUM)
├── Concrete07.h                       #
├── ConcreteCM.cpp                     # Confined/Unconfined (Priority: HIGH)
├── ConcreteCM.h                       #
├── ConcreteD.cpp                      # Damage model (Priority: MEDIUM)
├── ConcreteD.h                        #
├── ConcreteECThermal.cpp              # Eurocode thermal (Priority: LOW)
├── ConcreteECThermal.h                #
├── ConcretewBeta.cpp                  # With beta model (Priority: LOW)
├── ConcretewBeta.h                    #
├── FRPConfinedConcrete.cpp            # FRP confined (Priority: MEDIUM)
├── FRPConfinedConcrete.h              #
├── FRPConfinedConcrete02.cpp          # (Priority: MEDIUM)
└── FRPConfinedConcrete02.h            #

Concrete01 (Kent-Park-Scott):
- Parameters: fc, epsc0, fcu, epscu
- Parabolic ascending branch
- Linear descending branch
- No tensile strength

Concrete02 (with tension):
- Additional: ft, Ets
- Linear tension softening
- Tension stiffening
- Gap closing capability

Concrete04 (Popovics):
- n = Ec*epsc0/fc (curve parameter)
- Ascending: fc * n*(eps/epsc0) / (n-1 + (eps/epsc0)^n)
- Better fit to test data
- Unconfined and confined
```

#### K.2.3 Hysteretic Materials
```
SRC/material/uniaxial/
├── Hysteretic.cpp                     # General hysteretic (Priority: HIGH)
├── Hysteretic.h                       #
├── HystereticAsym.cpp                 # Asymmetric (Priority: MEDIUM)
├── HystereticAsym.h                   #
├── HystereticPoly.cpp                 # Polynomial backbone (Priority: LOW)
├── HystereticPoly.h                   #
├── HystereticSM.cpp                   # Sezen-Moehle (Priority: MEDIUM)
├── HystereticSM.h                     #
├── Bilin.cpp                          # Bilinear Ibarra (Priority: HIGH)
├── Bilin.h                            #
├── BilinMaterial.cpp                  # Modified bilinear (Priority: HIGH)
├── BilinMaterial.h                    #
├── IMKBilin.cpp                       # Ibarra-Medina-Krawinkler (Priority: HIGH)
├── IMKBilin.h                         #
├── IMKPeakOriented.cpp                # Peak-oriented (Priority: HIGH)
├── IMKPeakOriented.h                  #
├── IMKPinching.cpp                    # Pinching (Priority: HIGH)
├── IMKPinching.h                      #
├── ModIMKPeakOriented.cpp             # Modified IMK (Priority: MEDIUM)
├── ModIMKPeakOriented.h               #
├── ModIMKPeakOriented02.cpp           # (Priority: MEDIUM)
├── ModIMKPeakOriented02.h             #
├── ModIMKPinching.cpp                 # (Priority: MEDIUM)
├── ModIMKPinching.h                   #
├── ModIMKPinching02.cpp               # (Priority: MEDIUM)
├── ModIMKPinching02.h                 #
├── Pinching4.cpp                      # Pinching model (Priority: HIGH)
├── Pinching4.h                        #
├── PinchingLimitStateMaterial.cpp     # (Priority: MEDIUM)
└── PinchingLimitStateMaterial.h       #

Hysteretic material features:
- Multilinear backbone
- Strength/stiffness degradation
- Pinching behavior
- Energy-based damage index

IMK model parameters:
- Ke: elastic stiffness
- θy: yield rotation
- My: yield moment
- θcap: capping rotation
- θpc: post-capping rotation
- λS, λC, λK, λA: degradation parameters
- c: exponent for degradation
```

#### K.2.4 Damper Materials
```
SRC/material/uniaxial/
├── Viscous.cpp                        # Linear viscous (Priority: HIGH)
├── Viscous.h                          #
├── ViscousDamper.cpp                  # Nonlinear viscous (Priority: HIGH)
├── ViscousDamper.h                    #
├── Maxwell.cpp                        # Maxwell model (Priority: HIGH)
├── Maxwell.h                          #
├── KelvinVoigt.cpp                    # Kelvin-Voigt (Priority: HIGH)
├── KelvinVoigt.h                      #
├── DamperMaterial.cpp                 # General damper (Priority: MEDIUM)
├── DamperMaterial.h                   #
├── BilinearOilDamper.cpp              # Oil damper (Priority: MEDIUM)
├── BilinearOilDamper.h                #
├── SteelFractureDI.cpp                # Fracture damper (Priority: LOW)
├── SteelFractureDI.h                  #
├── BoucWenMaterial.cpp                # Bouc-Wen (Priority: HIGH)
├── BoucWenMaterial.h                  #
├── BoucWenOriginal.cpp                # Original formulation (Priority: MEDIUM)
├── BoucWenOriginal.h                  #
├── BWBN.cpp                           # Bouc-Wen-Baber-Noori (Priority: MEDIUM)
└── BWBN.h                             #

Viscous damper: F = C × |v|^α × sign(v)
- α = 1.0: linear
- α < 1.0: sublinear (typ. 0.3-0.5)
- α > 1.0: superlinear

Bouc-Wen model:
- z-dot = (A - (β + γ × sign(v × z)) × |z|^n) × v
- F = α × k × u + (1 - α) × k × z
- Parameters: α, ko, n, γ, β, Ao
```

#### K.2.5 nD (Multidimensional) Materials
```
SRC/material/nD/
├── ElasticIsotropicMaterial.cpp       # Elastic isotropic (Priority: EXISTS)
├── ElasticIsotropicMaterial.h         #
├── ElasticOrthotropicMaterial.cpp     # Elastic orthotropic (Priority: HIGH)
├── ElasticOrthotropicMaterial.h       #
├── J2Plasticity.cpp                   # J2 plasticity (Priority: EXISTS)
├── J2Plasticity.h                     #
├── J2PlaneStrain.cpp                  # Plane strain J2 (Priority: HIGH)
├── J2PlaneStrain.h                    #
├── J2PlaneStress.cpp                  # Plane stress J2 (Priority: HIGH)
├── J2PlaneStress.h                    #
├── J2Axisymm.cpp                      # Axisymmetric J2 (Priority: HIGH)
├── J2Axisymm.h                        #
├── J2CyclicBoundingSurface.cpp        # Bounding surface (Priority: MEDIUM)
├── J2CyclicBoundingSurface.h          #
├── DruckerPrager.cpp                  # Drucker-Prager (Priority: HIGH)
├── DruckerPrager.h                    #
├── DruckerPrager3D.cpp                # 3D DP (Priority: HIGH)
├── DruckerPrager3D.h                  #
├── DruckerPragerThermal.cpp           # With thermal (Priority: LOW)
├── DruckerPragerThermal.h             #
├── PressureIndependMultiYield.cpp     # Multi-yield surface (Priority: MEDIUM)
├── PressureIndependMultiYield.h       #
├── PressureDependMultiYield.cpp       # Pressure dependent (Priority: MEDIUM)
├── PressureDependMultiYield.h         #
├── PressureDependMultiYield02.cpp     # Updated PDMY (Priority: MEDIUM)
├── PressureDependMultiYield02.h       #
├── PressureDependMultiYield03.cpp     # (Priority: MEDIUM)
├── PressureDependMultiYield03.h       #
├── CycLiqCP.cpp                       # Cyclic liquefaction (Priority: LOW)
├── CycLiqCP.h                         #
├── CycLiqCPSP.cpp                     # (Priority: LOW)
├── CycLiqCPSP.h                       #
├── ManzariDafalias.cpp                # Manzari-Dafalias (Priority: MEDIUM)
├── ManzariDafalias.h                  #
├── ManzariDafalias3D.cpp              # 3D version (Priority: MEDIUM)
├── ManzariDafalias3D.h                #
├── ManzariDafaliasPlaneStrain.cpp     # Plane strain (Priority: MEDIUM)
├── ManzariDafaliasPlaneStrain.h       #
├── PM4Sand.cpp                        # PM4Sand model (Priority: MEDIUM)
├── PM4Sand.h                          #
├── PM4Silt.cpp                        # PM4Silt model (Priority: LOW)
├── PM4Silt.h                          #
├── AcousticMedium.cpp                 # Acoustic (Priority: LOW)
└── AcousticMedium.h                   #

Drucker-Prager yield function:
f = √J2 + α × I1 - k = 0

Where:
- J2 = second deviatoric stress invariant
- I1 = first stress invariant
- α = friction parameter
- k = cohesion parameter

Multi-yield surface concept:
- Nested yield surfaces for cyclic loading
- Kinematic hardening for each surface
- Captures modulus reduction with strain
```

#### K.2.6 Section Materials
```
SRC/material/section/
├── ElasticSection2d.cpp               # Elastic 2D section (Priority: HIGH)
├── ElasticSection2d.h                 #
├── ElasticSection3d.cpp               # Elastic 3D section (Priority: HIGH)
├── ElasticSection3d.h                 #
├── FiberSection2d.cpp                 # Fiber section 2D (Priority: CRITICAL)
├── FiberSection2d.h                   #
├── FiberSection3d.cpp                 # Fiber section 3D (Priority: HIGH)
├── FiberSection3d.h                   #
├── FiberSection2dThermal.cpp          # With thermal (Priority: LOW)
├── FiberSection2dThermal.h            #
├── FiberSection3dThermal.cpp          # (Priority: LOW)
├── FiberSection3dThermal.h            #
├── FiberSectionAsym3d.cpp             # Asymmetric 3D (Priority: MEDIUM)
├── FiberSectionAsym3d.h               #
├── FiberSectionWarping3d.cpp          # With warping (Priority: LOW)
├── FiberSectionWarping3d.h            #
├── MembranePlateFiberSection.cpp      # Shell fiber (Priority: MEDIUM)
├── MembranePlateFiberSection.h        #
├── LayeredShellFiberSection.cpp       # Layered shell (Priority: MEDIUM)
├── LayeredShellFiberSection.h         #
├── LayeredShellFiberSectionThermal.cpp # (Priority: LOW)
├── LayeredShellFiberSectionThermal.h   #
├── ElasticShearSection2d.cpp          # Shear section (Priority: MEDIUM)
├── ElasticShearSection2d.h            #
├── ElasticShearSection3d.cpp          # (Priority: MEDIUM)
├── ElasticShearSection3d.h            #
├── NDFiberSection2d.cpp               # nD fiber 2D (Priority: HIGH)
├── NDFiberSection2d.h                 #
├── NDFiberSection3d.cpp               # nD fiber 3D (Priority: HIGH)
├── NDFiberSection3d.h                 #
├── NDFiberSectionWarping2d.cpp        # With warping (Priority: LOW)
├── NDFiberSectionWarping2d.h          #
├── ParallelSection.cpp                # Parallel combination (Priority: LOW)
├── ParallelSection.h                  #
├── SectionAggregator.cpp              # Aggregated sections (Priority: HIGH)
├── SectionAggregator.h                #
├── WSection2d.cpp                     # Wide flange 2D (Priority: MEDIUM)
├── WSection2d.h                       #
├── WSection3d.cpp                     # Wide flange 3D (Priority: MEDIUM)
└── WSection3d.h                       #

Fiber section integration:
- Discretize cross-section into fibers
- Each fiber has: (y, z, A, material)
- Section response: Σ fibers
- State determination at fiber level
```

### K.3 Analysis Source Files

#### K.3.1 Integrator Source Files
```
SRC/analysis/integrator/
├── IncrementalIntegrator.cpp          # Base class (Priority: EXISTS)
├── IncrementalIntegrator.h            #
├── StaticIntegrator.cpp               # Static base (Priority: EXISTS)
├── StaticIntegrator.h                 #
├── TransientIntegrator.cpp            # Transient base (Priority: EXISTS)
├── TransientIntegrator.h              #
├── LoadControl.cpp                    # Load control (Priority: EXISTS)
├── LoadControl.h                      #
├── DisplacementControl.cpp            # Displacement control (Priority: HIGH)
├── DisplacementControl.h              #
├── ArcLength.cpp                      # Spherical arc-length (Priority: CRITICAL)
├── ArcLength.h                        #
├── ArcLength1.cpp                     # Crisfield method (Priority: HIGH)
├── ArcLength1.h                       #
├── MinUnbalDispNorm.cpp               # Minimum unbalanced (Priority: MEDIUM)
├── MinUnbalDispNorm.h                 #
├── DisplacementPath.cpp               # Displacement path (Priority: MEDIUM)
├── DisplacementPath.h                 #
├── HSConstraint.cpp                   # Hyper-spherical (Priority: LOW)
├── HSConstraint.h                     #
├── Newmark.cpp                        # Newmark-β (Priority: EXISTS)
├── Newmark.h                          #
├── Newmark1.cpp                       # Modified Newmark (Priority: MEDIUM)
├── Newmark1.h                         #
├── NewmarkExplicit.cpp                # Explicit Newmark (Priority: HIGH)
├── NewmarkExplicit.h                  #
├── NewmarkHSFixedNumIter.cpp          # Fixed iterations (Priority: LOW)
├── NewmarkHSFixedNumIter.h            #
├── NewmarkHSIncrDisp.cpp              # Incremental (Priority: LOW)
├── NewmarkHSIncrDisp.h                #
├── NewmarkHSIncrReduct.cpp            # With reduction (Priority: LOW)
├── NewmarkHSIncrReduct.h              #
├── HHT.cpp                            # Hilber-Hughes-Taylor (Priority: HIGH)
├── HHT.h                              #
├── HHT1.cpp                           # Modified HHT (Priority: MEDIUM)
├── HHT1.h                             #
├── HHTExplicit.cpp                    # Explicit HHT (Priority: MEDIUM)
├── HHTExplicit.h                      #
├── HHTGeneralized.cpp                 # Generalized HHT (Priority: LOW)
├── HHTGeneralized.h                   #
├── HHTGeneralizedExplicit.cpp         # (Priority: LOW)
├── HHTGeneralizedExplicit.h           #
├── HHTHSFixedNumIter.cpp              # (Priority: LOW)
├── HHTHSFixedNumIter.h                #
├── HHTHSIncrDisp.cpp                  # (Priority: LOW)
├── HHTHSIncrDisp.h                    #
├── HHTHSIncrReduct.cpp                # (Priority: LOW)
├── HHTHSIncrReduct.h                  #
├── WilsonTheta.cpp                    # Wilson-θ (Priority: MEDIUM)
├── WilsonTheta.h                      #
├── Collocation.cpp                    # Collocation method (Priority: MEDIUM)
├── Collocation.h                      #
├── CollocationHSFixedNumIter.cpp      # (Priority: LOW)
├── CollocationHSFixedNumIter.h        #
├── CollocationHSIncrDisp.cpp          # (Priority: LOW)
├── CollocationHSIncrDisp.h            #
├── CollocationHSIncrReduct.cpp        # (Priority: LOW)
├── CollocationHSIncrReduct.h          #
├── CentralDifference.cpp              # Central difference (Priority: HIGH)
├── CentralDifference.h                #
├── CentralDifferenceAlternative.cpp   # Alternative CD (Priority: LOW)
├── CentralDifferenceAlternative.h     #
├── CentralDifferenceNoDamping.cpp     # No damping CD (Priority: LOW)
├── CentralDifferenceNoDamping.h       #
├── AlphaOS.cpp                        # α-OS method (Priority: MEDIUM)
├── AlphaOS.h                          #
├── AlphaOSGeneralized.cpp             # Generalized α-OS (Priority: LOW)
├── AlphaOSGeneralized.h               #
├── GeneralizedAlpha.cpp               # Generalized-α (Priority: MEDIUM)
├── GeneralizedAlpha.h                 #
├── TRBDF2.cpp                         # TRBDF2 method (Priority: LOW)
├── TRBDF2.h                           #
├── Houbolt.cpp                        # Houbolt method (Priority: LOW)
├── Houbolt.h                          #
├── ParkLMS3.cpp                       # Park method (Priority: LOW)
├── ParkLMS3.h                         #
├── BackwardEuler.cpp                  # Backward Euler (Priority: MEDIUM)
├── BackwardEuler.h                    #
├── KRAlphaExplicit.cpp                # KR-α explicit (Priority: LOW)
├── KRAlphaExplicit.h                  #
├── StagedNewmark.cpp                  # Staged Newmark (Priority: LOW)
└── StagedNewmark.h                    #

Arc-length formulation:
- Constraint: ||Δu||² + ψ²Δλ² ||pref||² = Δs²
- Predictor: Δu = ±Δs × Kt⁻¹pref / ||Kt⁻¹pref||
- Corrector: solve augmented system
```

#### K.3.2 Algorithm Source Files
```
SRC/analysis/algorithm/
├── SolutionAlgorithm.cpp              # Base class (Priority: EXISTS)
├── SolutionAlgorithm.h                #
├── Linear.cpp                         # Linear solver (Priority: EXISTS)
├── Linear.h                           #
├── NewtonRaphson.cpp                  # Newton-Raphson (Priority: EXISTS)
├── NewtonRaphson.h                    #
├── ModifiedNewton.cpp                 # Modified Newton (Priority: HIGH)
├── ModifiedNewton.h                   #
├── NewtonLineSearch.cpp               # With line search (Priority: HIGH)
├── NewtonLineSearch.h                 #
├── Broyden.cpp                        # Broyden method (Priority: MEDIUM)
├── Broyden.h                          #
├── BFGS.cpp                           # BFGS method (Priority: MEDIUM)
├── BFGS.h                             #
├── KrylovNewton.cpp                   # Krylov-Newton (Priority: MEDIUM)
├── KrylovNewton.h                     #
├── PeriodicNewton.cpp                 # Periodic update (Priority: LOW)
├── PeriodicNewton.h                   #
├── AcceleratedNewton.cpp              # Accelerated (Priority: LOW)
├── AcceleratedNewton.h                #
├── RegulaFalsiLineSearch.cpp          # Regula falsi (Priority: MEDIUM)
├── RegulaFalsiLineSearch.h            #
├── BisectionLineSearch.cpp            # Bisection (Priority: MEDIUM)
├── BisectionLineSearch.h              #
├── SecantLineSearch.cpp               # Secant method (Priority: MEDIUM)
├── SecantLineSearch.h                 #
├── InitialInterpolatedLineSearch.cpp  # Interpolated (Priority: LOW)
└── InitialInterpolatedLineSearch.h    #

Line search algorithm:
- Find α such that g(α) = r(x + α×Δx)·Δx ≈ 0
- Armijo condition: f(x + α×Δx) ≤ f(x) + c₁×α×∇f·Δx
- Wolfe conditions add curvature check
```

#### K.3.3 System of Equations
```
SRC/system_of_eqn/
├── linearSOE/
│   ├── DomainSolver.cpp               # Domain decomposition (Priority: LOW)
│   ├── DomainSolver.h                 #
│   ├── bandGEN/                       # Banded general (Priority: EXISTS)
│   │   ├── BandGenLinLapackSolver.cpp
│   │   ├── BandGenLinLapackSolver.h
│   │   ├── BandGenLinSOE.cpp
│   │   └── BandGenLinSOE.h
│   ├── bandSPD/                       # Banded SPD (Priority: EXISTS)
│   │   ├── BandSPDLinLapackSolver.cpp
│   │   ├── BandSPDLinLapackSolver.h
│   │   ├── BandSPDLinSOE.cpp
│   │   └── BandSPDLinSOE.h
│   ├── fullGEN/                       # Full general (Priority: EXISTS)
│   │   ├── FullGenEigenSolver.cpp
│   │   ├── FullGenEigenSolver.h
│   │   ├── FullGenEigenSOE.cpp
│   │   └── FullGenEigenSOE.h
│   ├── profileSPD/                    # Profile SPD (Priority: HIGH)
│   │   ├── ProfileSPDLinDirectSkypackSolver.cpp
│   │   ├── ProfileSPDLinDirectSolver.cpp
│   │   ├── ProfileSPDLinSOE.cpp
│   │   └── ProfileSPDLinSOE.h
│   ├── sparseGEN/                     # Sparse general (Priority: HIGH)
│   │   ├── SparseGenColLinSOE.cpp
│   │   ├── SparseGenColLinSOE.h
│   │   ├── SuperLU.cpp
│   │   ├── SuperLU.h
│   │   ├── UmfpackGenLinSOE.cpp
│   │   └── UmfpackGenLinSolver.cpp
│   ├── sparseSYM/                     # Sparse symmetric (Priority: HIGH)
│   │   ├── SymSparseLinSOE.cpp
│   │   ├── SymSparseLinSOE.h
│   │   ├── SymSparseLinSolver.cpp
│   │   └── SymSparseLinSolver.h
│   ├── cg/                            # Conjugate gradient (Priority: HIGH)
│   │   ├── ConjugateGradientSolver.cpp
│   │   └── ConjugateGradientSolver.h
│   ├── diagonal/                      # Diagonal (Priority: HIGH)
│   │   ├── DiagonalDirectSolver.cpp
│   │   ├── DiagonalDirectSolver.h
│   │   ├── DiagonalSOE.cpp
│   │   └── DiagonalSOE.h
│   ├── mumps/                         # MUMPS solver (Priority: LOW)
│   │   ├── MumpsParallelSolver.cpp
│   │   ├── MumpsParallelSolver.h
│   │   ├── MumpsParallelSOE.cpp
│   │   └── MumpsParallelSOE.h
│   ├── petsc/                         # PETSc solver (Priority: LOW)
│   │   ├── PetscSolver.cpp
│   │   ├── PetscSolver.h
│   │   ├── PetscSOE.cpp
│   │   └── PetscSOE.h
│   └── itpack/                        # ITPACK solver (Priority: LOW)
│       ├── ItpackLinSolver.cpp
│       └── ItpackLinSOE.cpp
└── eigenSOE/                          # Eigenvalue problems (Priority: HIGH)
    ├── ArpackSOE.cpp                  # ARPACK interface
    ├── ArpackSOE.h
    ├── ArpackSolver.cpp
    ├── ArpackSolver.h
    ├── SymArpackSOE.cpp               # Symmetric ARPACK
    ├── SymArpackSOE.h
    ├── SymBandEigenSOE.cpp            # Band eigenvalue
    ├── SymBandEigenSOE.h
    ├── SymBandEigenSolver.cpp
    └── SymBandEigenSolver.h
```

---

## Appendix L: Complete CalculiX Element Catalog

### L.1 Three-Dimensional Solid Elements

```
Element Type     | Nodes | Integration | Notes
-----------------|-------|-------------|--------------------------------
C3D4            | 4     | 1 point     | Linear tetrahedron (EXISTS)
C3D4H           | 4     | 1 point     | Hybrid (incompressible)
C3D6            | 6     | 2 points    | Linear wedge/prism
C3D6H           | 6     | 2 points    | Hybrid wedge
C3D8            | 8     | 2×2×2       | Linear hexahedron (EXISTS)
C3D8H           | 8     | 2×2×2       | Hybrid hex
C3D8I           | 8     | 2×2×2       | Incompatible modes
C3D8R           | 8     | 1 point     | Reduced integration (Priority: HIGH)
C3D8RH          | 8     | 1 point     | Reduced hybrid
C3D10           | 10    | 4 points    | Quadratic tetrahedron (Priority: HIGH)
C3D10H          | 10    | 4 points    | Hybrid tet10
C3D10M          | 10    | 4 points    | Modified tet10
C3D10MH         | 10    | 4 points    | Modified hybrid
C3D15           | 15    | 9 points    | Quadratic wedge
C3D15H          | 15    | 9 points    | Hybrid wedge15
C3D20           | 20    | 3×3×3       | Quadratic hexahedron (Priority: HIGH)
C3D20H          | 20    | 3×3×3       | Hybrid hex20
C3D20R          | 20    | 2×2×2       | Reduced hex20 (Priority: HIGH)
C3D20RH         | 20    | 2×2×2       | Reduced hybrid
C3D27           | 27    | 3×3×3       | Lagrangian hex
C3D27H          | 27    | 3×3×3       | Hybrid hex27
C3D27R          | 27    | 2×2×2       | Reduced hex27
```

### L.2 Plane Stress Elements

```
Element Type     | Nodes | Integration | Notes
-----------------|-------|-------------|--------------------------------
CPS3            | 3     | 1 point     | Linear triangle (EXISTS)
CPS4            | 4     | 2×2         | Linear quadrilateral (EXISTS)
CPS4I           | 4     | 2×2         | Incompatible modes (Priority: HIGH)
CPS4R           | 4     | 1 point     | Reduced integration (Priority: HIGH)
CPS6            | 6     | 3 points    | Quadratic triangle (Priority: HIGH)
CPS6M           | 6     | 3 points    | Modified
CPS8            | 8     | 3×3         | Quadratic quad (Priority: HIGH)
CPS8R           | 8     | 2×2         | Reduced quad8 (Priority: HIGH)
```

### L.3 Plane Strain Elements

```
Element Type     | Nodes | Integration | Notes
-----------------|-------|-------------|--------------------------------
CPE3            | 3     | 1 point     | Linear triangle
CPE3H           | 3     | 1 point     | Hybrid
CPE4            | 4     | 2×2         | Linear quad
CPE4H           | 4     | 2×2         | Hybrid
CPE4I           | 4     | 2×2         | Incompatible modes
CPE4R           | 4     | 1 point     | Reduced (Priority: HIGH)
CPE4RH          | 4     | 1 point     | Reduced hybrid
CPE6            | 6     | 3 points    | Quadratic triangle (Priority: HIGH)
CPE6H           | 6     | 3 points    | Hybrid
CPE6M           | 6     | 3 points    | Modified
CPE6MH          | 6     | 3 points    | Modified hybrid
CPE8            | 8     | 3×3         | Quadratic quad (Priority: HIGH)
CPE8H           | 8     | 3×3         | Hybrid
CPE8R           | 8     | 2×2         | Reduced (Priority: HIGH)
CPE8RH          | 8     | 2×2         | Reduced hybrid
```

### L.4 Axisymmetric Elements

```
Element Type     | Nodes | Integration | Notes
-----------------|-------|-------------|--------------------------------
CAX3            | 3     | 1 point     | Linear triangle (Priority: HIGH)
CAX3H           | 3     | 1 point     | Hybrid
CAX4            | 4     | 2×2         | Linear quad (Priority: HIGH)
CAX4H           | 4     | 2×2         | Hybrid
CAX4I           | 4     | 2×2         | Incompatible modes
CAX4R           | 4     | 1 point     | Reduced (Priority: HIGH)
CAX4RH          | 4     | 1 point     | Reduced hybrid
CAX6            | 6     | 3 points    | Quadratic triangle (Priority: MEDIUM)
CAX6H           | 6     | 3 points    | Hybrid
CAX6M           | 6     | 3 points    | Modified
CAX6MH          | 6     | 3 points    | Modified hybrid
CAX8            | 8     | 3×3         | Quadratic quad (Priority: MEDIUM)
CAX8H           | 8     | 3×3         | Hybrid
CAX8R           | 8     | 2×2         | Reduced (Priority: MEDIUM)
CAX8RH          | 8     | 2×2         | Reduced hybrid
```

### L.5 Shell Elements

```
Element Type     | Nodes | Integration | Notes
-----------------|-------|-------------|--------------------------------
S3              | 3     | 1 in-plane  | Linear triangle (Priority: HIGH)
S3R             | 3     | 1           | Reduced
S4              | 4     | 2×2         | Linear quad (Priority: HIGH)
S4R             | 4     | 1           | Reduced (Priority: HIGH)
S4R5            | 4     | 1           | 5 DOF per node
S6              | 6     | 3 points    | Quadratic triangle (Priority: MEDIUM)
S6R             | 6     | 1           | Reduced
S8              | 8     | 3×3         | Quadratic quad (Priority: MEDIUM)
S8R             | 8     | 2×2         | Reduced (Priority: MEDIUM)
S8R5            | 8     | 2×2         | 5 DOF per node
STRI3           | 3     | 1           | Thin shell triangle
STRI65          | 6     | 3           | Thin shell 6-node
```

### L.6 Beam Elements

```
Element Type     | Nodes | Integration | Notes
-----------------|-------|-------------|--------------------------------
B21             | 2     | 1 point     | 2D linear beam (Priority: HIGH)
B22             | 3     | 2 points    | 2D quadratic (Priority: MEDIUM)
B23             | 2     | 2 points    | 2D cubic Hermite (Priority: HIGH)
B31             | 2     | 1 point     | 3D linear beam (Priority: HIGH)
B31OS           | 2     | 1 point     | Open section
B31H            | 2     | 1 point     | Hybrid
B32             | 3     | 2 points    | 3D quadratic (Priority: MEDIUM)
B32OS           | 3     | 2 points    | Open section
B32H            | 3     | 2 points    | Hybrid
B33             | 2     | 2 points    | 3D cubic (Priority: HIGH)
B33H            | 2     | 2 points    | Hybrid
```

### L.7 Special Elements

```
Element Type     | Description                    | Priority
-----------------|--------------------------------|----------
CONN2D2         | 2D connector                   | HIGH
CONN3D2         | 3D connector                   | HIGH
DCOUP2D         | 2D distributed coupling        | MEDIUM
DCOUP3D         | 3D distributed coupling        | MEDIUM
GAPUNI          | Gap element                    | HIGH
GAPCYL          | Cylindrical gap                | MEDIUM
GAPSPH          | Spherical gap                  | MEDIUM
DASHPOTA        | Axial dashpot                  | HIGH
DASHPOT1        | 1-node dashpot                 | HIGH
DASHPOT2        | 2-node dashpot                 | HIGH
SPRING1         | 1-node spring                  | HIGH
SPRING2         | 2-node spring                  | HIGH
SPRINGA         | Axial spring                   | HIGH
MASS            | Point mass                     | HIGH
ROTARYI         | Rotary inertia                 | HIGH
```

---

## Appendix M: Detailed API Specifications

### M.1 Element API Requirements

Each element must implement the following interface:

```python
class ElementBase(ABC):
    """
    Abstract base class for finite elements.
    
    All femlabpy elements must inherit from this class and implement
    the required abstract methods.
    """
    
    # Class attributes
    ndim: int           # Spatial dimension (1, 2, or 3)
    nnodes: int         # Number of nodes
    ndof_per_node: int  # DOFs per node
    ndof: int           # Total DOFs = nnodes × ndof_per_node
    nip: int            # Number of integration points
    
    @abstractmethod
    def __init__(self, nodes: np.ndarray, properties: dict):
        """
        Initialize element.
        
        Parameters
        ----------
        nodes : ndarray, shape (nnodes, ndim)
            Node coordinates in global system
        properties : dict
            Element properties (material, geometry, etc.)
        """
        pass
    
    @abstractmethod
    def stiffness_matrix(self) -> np.ndarray:
        """
        Compute element stiffness matrix in global coordinates.
        
        Returns
        -------
        K : ndarray, shape (ndof, ndof)
            Element stiffness matrix
        """
        pass
    
    @abstractmethod
    def mass_matrix(self, lumped: bool = False) -> np.ndarray:
        """
        Compute element mass matrix in global coordinates.
        
        Parameters
        ----------
        lumped : bool, optional
            If True, return lumped (diagonal) mass matrix.
            Default is False (consistent mass).
        
        Returns
        -------
        M : ndarray, shape (ndof, ndof)
            Element mass matrix
        """
        pass
    
    @abstractmethod
    def internal_force(self, u: np.ndarray) -> np.ndarray:
        """
        Compute element internal force vector.
        
        Parameters
        ----------
        u : ndarray, shape (ndof,)
            Element displacement vector in global coordinates
        
        Returns
        -------
        f : ndarray, shape (ndof,)
            Element internal force vector
        """
        pass
    
    @abstractmethod
    def strain_at_points(self, u: np.ndarray) -> np.ndarray:
        """
        Compute strains at integration points.
        
        Parameters
        ----------
        u : ndarray, shape (ndof,)
            Element displacement vector
        
        Returns
        -------
        epsilon : ndarray, shape (nip, nstrain)
            Strains at integration points
            nstrain depends on element type:
            - 1D: 1 (axial strain)
            - 2D: 3 (εxx, εyy, γxy)
            - 3D: 6 (εxx, εyy, εzz, γxy, γyz, γxz)
        """
        pass
    
    @abstractmethod
    def stress_at_points(self, u: np.ndarray) -> np.ndarray:
        """
        Compute stresses at integration points.
        
        Parameters
        ----------
        u : ndarray, shape (ndof,)
            Element displacement vector
        
        Returns
        -------
        sigma : ndarray, shape (nip, nstress)
            Stresses at integration points
        """
        pass
    
    @abstractmethod
    def transformation_matrix(self) -> np.ndarray:
        """
        Compute transformation matrix from local to global coordinates.
        
        Returns
        -------
        T : ndarray, shape (ndof, ndof)
            Transformation matrix
        """
        pass
    
    def geometric_stiffness(self, u: np.ndarray) -> np.ndarray:
        """
        Compute geometric stiffness matrix for nonlinear analysis.
        
        Parameters
        ----------
        u : ndarray, shape (ndof,)
            Element displacement vector
        
        Returns
        -------
        Kg : ndarray, shape (ndof, ndof)
            Geometric stiffness matrix
            
        Note: Default implementation returns zero matrix.
        Override for elements supporting geometric nonlinearity.
        """
        return np.zeros((self.ndof, self.ndof))
    
    def damping_matrix(self, alpha: float, beta: float) -> np.ndarray:
        """
        Compute element damping matrix using Rayleigh damping.
        
        C = α × M + β × K
        
        Parameters
        ----------
        alpha : float
            Mass proportional damping coefficient
        beta : float
            Stiffness proportional damping coefficient
        
        Returns
        -------
        C : ndarray, shape (ndof, ndof)
            Element damping matrix
        """
        return alpha * self.mass_matrix() + beta * self.stiffness_matrix()
    
    def dof_indices(self) -> np.ndarray:
        """
        Return global DOF indices for element assembly.
        
        Returns
        -------
        indices : ndarray, shape (ndof,)
            Global DOF indices
        """
        # Implementation depends on mesh connectivity
        pass
```

### M.2 Material API Requirements

```python
class UniaxialMaterialBase(ABC):
    """
    Abstract base class for uniaxial materials.
    """
    
    @abstractmethod
    def __init__(self, **kwargs):
        """Initialize material with parameters."""
        pass
    
    @abstractmethod
    def get_stress(self, strain: float) -> Tuple[float, float]:
        """
        Compute stress and tangent modulus for given strain.
        
        Parameters
        ----------
        strain : float
            Current strain value
        
        Returns
        -------
        stress : float
            Current stress
        Et : float
            Tangent modulus (dσ/dε)
        """
        pass
    
    @abstractmethod
    def commit_state(self) -> None:
        """
        Commit the current state as converged.
        Called after successful equilibrium iteration.
        """
        pass
    
    @abstractmethod
    def revert_to_last_commit(self) -> None:
        """
        Revert to the last committed state.
        Called when iteration fails to converge.
        """
        pass
    
    @abstractmethod
    def revert_to_start(self) -> None:
        """
        Revert to initial (undeformed) state.
        """
        pass
    
    @abstractmethod
    def copy(self) -> 'UniaxialMaterialBase':
        """
        Create a copy of the material with independent state.
        
        Returns
        -------
        material : UniaxialMaterialBase
            Copy of material
        """
        pass
    
    @property
    @abstractmethod
    def strain(self) -> float:
        """Current strain value."""
        pass
    
    @property
    @abstractmethod
    def stress(self) -> float:
        """Current stress value."""
        pass
    
    @property
    @abstractmethod
    def tangent(self) -> float:
        """Current tangent modulus."""
        pass


class NDMaterialBase(ABC):
    """
    Abstract base class for multidimensional materials.
    """
    
    # Class attributes
    ndim: int           # Spatial dimension (2 or 3)
    nstrain: int        # Number of strain components
    nstress: int        # Number of stress components
    
    @abstractmethod
    def __init__(self, **kwargs):
        """Initialize material with parameters."""
        pass
    
    @abstractmethod
    def get_stress(self, strain: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute stress and tangent stiffness for given strain.
        
        Parameters
        ----------
        strain : ndarray, shape (nstrain,)
            Current strain vector
            2D: [εxx, εyy, γxy]
            3D: [εxx, εyy, εzz, γxy, γyz, γxz]
        
        Returns
        -------
        stress : ndarray, shape (nstress,)
            Current stress vector
        Ct : ndarray, shape (nstress, nstrain)
            Tangent constitutive matrix
        """
        pass
    
    @abstractmethod
    def commit_state(self) -> None:
        """Commit the current state as converged."""
        pass
    
    @abstractmethod
    def revert_to_last_commit(self) -> None:
        """Revert to the last committed state."""
        pass
    
    @abstractmethod
    def revert_to_start(self) -> None:
        """Revert to initial state."""
        pass
    
    @abstractmethod
    def copy(self) -> 'NDMaterialBase':
        """Create a copy of the material."""
        pass


class SectionBase(ABC):
    """
    Abstract base class for beam/shell sections.
    """
    
    @abstractmethod
    def __init__(self, **kwargs):
        """Initialize section."""
        pass
    
    @abstractmethod
    def get_section_response(self, phi: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute section forces and stiffness for given deformations.
        
        Parameters
        ----------
        phi : ndarray
            Section deformations
            2D beam: [ε₀, κ] (axial strain, curvature)
            3D beam: [ε₀, κy, κz] or [ε₀, κy, κz, γy, γz, θ]
        
        Returns
        -------
        s : ndarray
            Section forces [N, M] or [N, My, Mz, ...]
        ks : ndarray
            Section stiffness matrix
        """
        pass
    
    @abstractmethod
    def commit_state(self) -> None:
        """Commit the current state."""
        pass
    
    @abstractmethod
    def revert_to_last_commit(self) -> None:
        """Revert to last committed state."""
        pass
```

### M.3 Solver API Requirements

```python
class StaticSolverBase(ABC):
    """Base class for static solvers."""
    
    @abstractmethod
    def solve(
        self,
        K: sparse.csr_matrix,
        f: np.ndarray,
        bc: dict
    ) -> np.ndarray:
        """
        Solve the static system K × u = f with boundary conditions.
        
        Parameters
        ----------
        K : sparse matrix, shape (ndof, ndof)
            Global stiffness matrix
        f : ndarray, shape (ndof,)
            Global force vector
        bc : dict
            Boundary conditions {dof_index: prescribed_value}
        
        Returns
        -------
        u : ndarray, shape (ndof,)
            Displacement solution
        """
        pass


class NonlinearSolverBase(ABC):
    """Base class for nonlinear solvers."""
    
    @abstractmethod
    def solve(
        self,
        model,
        load_pattern: Callable,
        max_steps: int = 100,
        tolerance: float = 1e-6
    ) -> dict:
        """
        Solve nonlinear static problem.
        
        Parameters
        ----------
        model : FiniteElementModel
            FE model with elements and BCs
        load_pattern : callable
            Function returning load vector for given load factor
        max_steps : int
            Maximum number of load steps
        tolerance : float
            Convergence tolerance
        
        Returns
        -------
        results : dict
            {
                'displacements': ndarray,
                'reactions': ndarray,
                'load_factors': ndarray,
                'converged': bool
            }
        """
        pass


class DynamicSolverBase(ABC):
    """Base class for dynamic solvers."""
    
    @abstractmethod
    def solve(
        self,
        M: sparse.csr_matrix,
        C: sparse.csr_matrix,
        K: sparse.csr_matrix,
        f: Callable,
        u0: np.ndarray,
        v0: np.ndarray,
        dt: float,
        nsteps: int,
        bc: dict
    ) -> dict:
        """
        Solve the dynamic system M×ü + C×u̇ + K×u = f(t).
        
        Parameters
        ----------
        M : sparse matrix
            Mass matrix
        C : sparse matrix
            Damping matrix
        K : sparse matrix
            Stiffness matrix
        f : callable
            Time-dependent force function f(t) -> ndarray
        u0 : ndarray
            Initial displacement
        v0 : ndarray
            Initial velocity
        dt : float
            Time step
        nsteps : int
            Number of time steps
        bc : dict
            Boundary conditions
        
        Returns
        -------
        results : dict
            {
                'time': ndarray,
                'displacement': ndarray,
                'velocity': ndarray,
                'acceleration': ndarray
            }
        """
        pass


class ModalSolverBase(ABC):
    """Base class for modal solvers."""
    
    @abstractmethod
    def solve(
        self,
        K: sparse.csr_matrix,
        M: sparse.csr_matrix,
        n_modes: int = 10,
        bc: dict = None
    ) -> dict:
        """
        Solve the eigenvalue problem K×φ = ω²×M×φ.
        
        Parameters
        ----------
        K : sparse matrix
            Stiffness matrix
        M : sparse matrix
            Mass matrix
        n_modes : int
            Number of modes to compute
        bc : dict, optional
            Boundary conditions (constrained DOFs)
        
        Returns
        -------
        results : dict
            {
                'eigenvalues': ndarray, shape (n_modes,)
                'frequencies': ndarray, shape (n_modes,)  # Hz
                'periods': ndarray, shape (n_modes,)      # seconds
                'mode_shapes': ndarray, shape (ndof, n_modes)
            }
        """
        pass
```

---

## Appendix N: Performance Benchmarks Reference

### N.1 Assembly Performance Targets

| Operation | Current | Target | Notes |
|-----------|---------|--------|-------|
| T3 element K_e | 0.1 ms | 0.05 ms | Vectorized |
| Q4 element K_e | 0.2 ms | 0.1 ms | Vectorized |
| H8 element K_e | 0.5 ms | 0.2 ms | Vectorized |
| T10 element K_e | N/A | 0.3 ms | New |
| H20 element K_e | N/A | 0.8 ms | New |
| Global assembly (10K) | 2.0 s | 0.5 s | Sparse COO |
| Global assembly (100K) | 25 s | 5 s | Sparse COO |

### N.2 Solver Performance Targets

| Solver | Problem Size | Current | Target |
|--------|-------------|---------|--------|
| Direct (UMFPACK) | 10K DOF | 0.3 s | 0.2 s |
| Direct (UMFPACK) | 100K DOF | 15 s | 10 s |
| Iterative (CG) | 10K DOF | 0.5 s | 0.3 s |
| Iterative (CG) | 100K DOF | 8 s | 5 s |
| Iterative (CG) | 1M DOF | N/A | 60 s |
| Modal (ARPACK) | 10 modes, 10K | 0.5 s | 0.3 s |
| Modal (ARPACK) | 50 modes, 100K | 30 s | 15 s |

### N.3 Memory Usage Targets

| Data Structure | Size | Current | Target |
|----------------|------|---------|--------|
| Sparse K (10K) | N/A | 50 MB | 20 MB |
| Sparse K (100K) | N/A | 800 MB | 300 MB |
| Element storage | per elem | 1 KB | 0.5 KB |
| Material state | per IP | 200 B | 100 B |

---

## Appendix O: Code Organization Standards

### O.1 Module Structure

```
femlabpy/
├── __init__.py                 # Public API exports
├── _helpers.py                 # Internal utilities
├── _version.py                 # Version info
├── types.py                    # Type definitions
│
├── elements/                   # Element implementations
│   ├── __init__.py             # Element exports
│   ├── _base.py                # ElementBase class
│   ├── _gauss.py               # Gaussian quadrature
│   ├── _shape.py               # Shape function library
│   ├── triangles.py            # T3, T6
│   ├── quads.py                # Q4, Q8, Q9
│   ├── bars.py                 # Bar2, Bar3
│   ├── solids.py               # T4, T10, H8, H20
│   ├── beams2d.py              # Beam2d (NEW)
│   ├── beams3d.py              # Beam3d (NEW)
│   ├── shells.py               # DKT, MITC4 (NEW)
│   └── springs.py              # ZeroLength (NEW)
│
├── materials/                  # Material models
│   ├── __init__.py             # Material exports
│   ├── _base.py                # MaterialBase classes
│   ├── elastic.py              # Elastic materials
│   ├── plasticity.py           # J2 plasticity (EXISTS)
│   ├── invariants.py           # Stress invariants
│   ├── steel.py                # Steel01, Steel02 (NEW)
│   ├── concrete.py             # Concrete01, Concrete02 (NEW)
│   ├── hysteretic.py           # Hysteretic, Pinching4 (NEW)
│   ├── dampers.py              # Viscous, Maxwell (NEW)
│   └── soil.py                 # DruckerPrager, PDMY (NEW)
│
├── sections/                   # Section models (NEW)
│   ├── __init__.py
│   ├── _base.py                # SectionBase
│   ├── elastic.py              # ElasticSection
│   ├── fiber.py                # FiberSection
│   └── shapes.py               # Section shape definitions
│
├── analysis/                   # Analysis types (NEW)
│   ├── __init__.py
│   ├── linear.py               # Linear static (EXISTS)
│   ├── nonlinear.py            # Nonlinear static (NEW)
│   ├── arclength.py            # Arc-length (NEW)
│   ├── buckling.py             # Linear buckling (NEW)
│   ├── spectrum.py             # Response spectrum (NEW)
│   └── thermal.py              # Thermal analysis (NEW)
│
├── assembly.py                 # Global assembly
├── boundary.py                 # Boundary conditions
├── loads.py                    # Load patterns
├── solvers.py                  # Linear solvers
├── modal.py                    # Eigenvalue analysis
├── dynamics.py                 # Time history
├── damping.py                  # Damping models
├── periodic.py                 # Periodic analysis
│
├── io/                         # Input/Output (NEW)
│   ├── __init__.py
│   ├── gmsh.py                 # Gmsh interface
│   ├── opensees.py             # OpenSees format
│   ├── calculix.py             # CalculiX format
│   ├── vtk.py                  # VTK export
│   └── hdf5.py                 # HDF5 storage
│
└── utils/                      # Utilities (NEW)
    ├── __init__.py
    ├── mesh.py                 # Mesh utilities
    ├── postprocess.py          # Post-processing
    └── visualization.py        # Plotting
```

### O.2 Naming Conventions

| Item | Convention | Example |
|------|------------|---------|
| Modules | lowercase_underscore | `beam_elements.py` |
| Classes | PascalCase | `ElasticBeam2d` |
| Functions | lowercase_underscore | `compute_stiffness` |
| Constants | UPPERCASE_UNDERSCORE | `MAX_ITERATIONS` |
| Private | _leading_underscore | `_internal_helper` |
| Type vars | SingleCapital | `T`, `MaterialT` |

### O.3 Documentation Standards

```python
def compute_element_stiffness(
    xe: np.ndarray,
    ge: np.ndarray,
    element_type: str = "Q4"
) -> np.ndarray:
    """
    Compute element stiffness matrix.
    
    Mathematical Formulation
    ------------------------
    The stiffness matrix is computed as:
    
    .. math::
        K_e = \\int_\\Omega B^T D B \\, d\\Omega
    
    where B is the strain-displacement matrix and D is the
    constitutive matrix.
    
    Parameters
    ----------
    xe : ndarray, shape (nnodes, ndim)
        Node coordinates for the element.
    ge : ndarray
        Material/geometric properties.
        For plane stress: [E, nu, t] (Young's modulus, Poisson's ratio, thickness)
        For plane strain: [E, nu]
    element_type : str, optional
        Element type identifier. Default is "Q4".
        Supported: "T3", "T6", "Q4", "Q8", "Q9"
    
    Returns
    -------
    ke : ndarray, shape (ndof, ndof)
        Element stiffness matrix in global coordinates.
    
    Raises
    ------
    ValueError
        If element_type is not supported.
    
    Examples
    --------
    >>> xe = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    >>> ge = np.array([210e9, 0.3, 0.01])  # E, nu, t
    >>> ke = compute_element_stiffness(xe, ge, "Q4")
    >>> ke.shape
    (8, 8)
    
    See Also
    --------
    compute_mass_matrix : Compute element mass matrix.
    assemble_global : Assemble global matrices.
    
    References
    ----------
    .. [1] Bathe, K.J. "Finite Element Procedures", 2nd ed., 2014.
    .. [2] Zienkiewicz, O.C. and Taylor, R.L. "The Finite Element Method",
           7th ed., 2013.
    """
    pass
```

---

## Appendix P: Numerical Integration Rules

### P.1 One-Dimensional Gauss-Legendre Quadrature

```
n=1:  ξ₁ = 0.0                          w₁ = 2.0

n=2:  ξ₁ = -1/√3 ≈ -0.577350269        w₁ = 1.0
      ξ₂ = +1/√3 ≈ +0.577350269        w₂ = 1.0

n=3:  ξ₁ = -√(3/5) ≈ -0.774596669      w₁ = 5/9 ≈ 0.555555556
      ξ₂ = 0.0                          w₂ = 8/9 ≈ 0.888888889
      ξ₃ = +√(3/5) ≈ +0.774596669      w₃ = 5/9 ≈ 0.555555556

n=4:  ξ₁ = -0.861136311594953          w₁ = 0.347854845137454
      ξ₂ = -0.339981043584856          w₂ = 0.652145154862546
      ξ₃ = +0.339981043584856          w₃ = 0.652145154862546
      ξ₄ = +0.861136311594953          w₄ = 0.347854845137454

n=5:  ξ₁ = -0.906179845938664          w₁ = 0.236926885056189
      ξ₂ = -0.538469310105683          w₂ = 0.478628670499366
      ξ₃ = 0.0                          w₃ = 0.568888888888889
      ξ₄ = +0.538469310105683          w₄ = 0.478628670499366
      ξ₅ = +0.906179845938664          w₅ = 0.236926885056189
```

### P.2 Triangular Integration Rules

#### 1-Point Rule (Order 1)
```
Point 1: (1/3, 1/3, 1/3)  weight = 1.0
```

#### 3-Point Rule (Order 2)
```
Point 1: (2/3, 1/6, 1/6)  weight = 1/3
Point 2: (1/6, 2/3, 1/6)  weight = 1/3
Point 3: (1/6, 1/6, 2/3)  weight = 1/3
```

#### 4-Point Rule (Order 3)
```
Point 1: (1/3, 1/3, 1/3)  weight = -27/48
Point 2: (3/5, 1/5, 1/5)  weight = 25/48
Point 3: (1/5, 3/5, 1/5)  weight = 25/48
Point 4: (1/5, 1/5, 3/5)  weight = 25/48
```

#### 6-Point Rule (Order 4)
```
α = 0.816847572980459
β = 0.091576213509771
γ = 0.108103018168070
δ = 0.445948490915965

Point 1: (α, β, β)  weight = 0.109951743655322
Point 2: (β, α, β)  weight = 0.109951743655322
Point 3: (β, β, α)  weight = 0.109951743655322
Point 4: (γ, δ, δ)  weight = 0.223381589678011
Point 5: (δ, γ, δ)  weight = 0.223381589678011
Point 6: (δ, δ, γ)  weight = 0.223381589678011
```

#### 7-Point Rule (Order 5)
```
α = 0.797426985353087
β = 0.101286507323456
γ = 0.470142064105115
δ = 0.059715871789770

Point 1: (1/3, 1/3, 1/3)    weight = 0.225000000000000
Point 2: (α, β, β)          weight = 0.125939180544827
Point 3: (β, α, β)          weight = 0.125939180544827
Point 4: (β, β, α)          weight = 0.125939180544827
Point 5: (γ, γ, δ)          weight = 0.132394152788506
Point 6: (γ, δ, γ)          weight = 0.132394152788506
Point 7: (δ, γ, γ)          weight = 0.132394152788506
```

### P.3 Tetrahedral Integration Rules

#### 1-Point Rule (Order 1)
```
Point 1: (1/4, 1/4, 1/4, 1/4)  weight = 1.0
```

#### 4-Point Rule (Order 2)
```
α = (5 - √5) / 20 ≈ 0.138196601125011
β = (5 + 3√5) / 20 ≈ 0.585410196624968

Point 1: (α, α, α, β)  weight = 1/4
Point 2: (α, α, β, α)  weight = 1/4
Point 3: (α, β, α, α)  weight = 1/4
Point 4: (β, α, α, α)  weight = 1/4
```

#### 5-Point Rule (Order 3)
```
Point 1: (1/4, 1/4, 1/4, 1/4)  weight = -4/5
Point 2: (1/2, 1/6, 1/6, 1/6)  weight = 9/20
Point 3: (1/6, 1/2, 1/6, 1/6)  weight = 9/20
Point 4: (1/6, 1/6, 1/2, 1/6)  weight = 9/20
Point 5: (1/6, 1/6, 1/6, 1/2)  weight = 9/20
```

#### 11-Point Rule (Order 4)
```
α = (1 + √(5/14)) / 4 ≈ 0.399403576166799
β = (1 - √(5/14)) / 4 ≈ 0.100596423833201

Point 1:  (1/4, 1/4, 1/4, 1/4)  weight = -74/5625
Points 2-5:   (α, β, β, β) and permutations, weight = 343/45000
Points 6-9:   (β, α, α, α) and permutations, weight = 343/45000
Points 10-11: (1/2, 1/6, 1/6, 1/6) + edge midpoints, weight = 56/2250
```

### P.4 Gauss-Lobatto Quadrature (for beam elements)

```
n=2:  ξ₁ = -1.0                         w₁ = 1.0
      ξ₂ = +1.0                         w₂ = 1.0

n=3:  ξ₁ = -1.0                         w₁ = 1/3
      ξ₂ = 0.0                          w₂ = 4/3
      ξ₃ = +1.0                         w₃ = 1/3

n=4:  ξ₁ = -1.0                         w₁ = 1/6
      ξ₂ = -√(1/5) ≈ -0.447213595       w₂ = 5/6
      ξ₃ = +√(1/5) ≈ +0.447213595       w₃ = 5/6
      ξ₄ = +1.0                         w₄ = 1/6

n=5:  ξ₁ = -1.0                         w₁ = 1/10
      ξ₂ = -√(3/7) ≈ -0.654653671       w₂ = 49/90
      ξ₃ = 0.0                          w₃ = 32/45
      ξ₄ = +√(3/7) ≈ +0.654653671       w₄ = 49/90
      ξ₅ = +1.0                         w₅ = 1/10
```

### P.5 Implementation in Python

```python
def gauss_legendre_1d(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Gauss-Legendre quadrature points and weights.
    
    Parameters
    ----------
    n : int
        Number of integration points (1 to 10 supported)
    
    Returns
    -------
    xi : ndarray, shape (n,)
        Integration points in [-1, 1]
    w : ndarray, shape (n,)
        Integration weights
    
    Examples
    --------
    >>> xi, w = gauss_legendre_1d(3)
    >>> xi
    array([-0.77459667,  0.        ,  0.77459667])
    >>> w
    array([0.55555556, 0.88888889, 0.55555556])
    """
    if n == 1:
        return np.array([0.0]), np.array([2.0])
    elif n == 2:
        s = 1.0 / np.sqrt(3.0)
        return np.array([-s, s]), np.array([1.0, 1.0])
    elif n == 3:
        s = np.sqrt(3.0/5.0)
        return np.array([-s, 0.0, s]), np.array([5/9, 8/9, 5/9])
    # ... continue for n=4 through n=10
    else:
        # Use numpy's polynomial roots
        return np.polynomial.legendre.leggauss(n)


def gauss_triangle(order: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Gauss quadrature for triangles in area coordinates.
    
    Parameters
    ----------
    order : int
        Polynomial order (1 to 5)
    
    Returns
    -------
    L : ndarray, shape (nip, 3)
        Area coordinates [L1, L2, L3] for each integration point
    w : ndarray, shape (nip,)
        Integration weights (sum to 1.0 for unit triangle)
    """
    if order == 1:
        # 1-point rule
        L = np.array([[1/3, 1/3, 1/3]])
        w = np.array([1.0])
    elif order == 2:
        # 3-point rule
        L = np.array([
            [2/3, 1/6, 1/6],
            [1/6, 2/3, 1/6],
            [1/6, 1/6, 2/3]
        ])
        w = np.array([1/3, 1/3, 1/3])
    elif order == 3:
        # 4-point rule
        L = np.array([
            [1/3, 1/3, 1/3],
            [3/5, 1/5, 1/5],
            [1/5, 3/5, 1/5],
            [1/5, 1/5, 3/5]
        ])
        w = np.array([-27/48, 25/48, 25/48, 25/48])
    # ... continue for higher orders
    
    return L, w


def gauss_tetrahedron(order: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Gauss quadrature for tetrahedra in volume coordinates.
    
    Parameters
    ----------
    order : int
        Polynomial order (1 to 4)
    
    Returns
    -------
    L : ndarray, shape (nip, 4)
        Volume coordinates [L1, L2, L3, L4] for each point
    w : ndarray, shape (nip,)
        Integration weights (sum to 1.0 for unit tetrahedron)
    """
    if order == 1:
        # 1-point rule
        L = np.array([[0.25, 0.25, 0.25, 0.25]])
        w = np.array([1.0])
    elif order == 2:
        # 4-point rule
        a = 0.138196601125011
        b = 0.585410196624968
        L = np.array([
            [a, a, a, b],
            [a, a, b, a],
            [a, b, a, a],
            [b, a, a, a]
        ])
        w = np.array([0.25, 0.25, 0.25, 0.25])
    # ... continue for higher orders
    
    return L, w
```

---

## Appendix Q: Shape Function Library

### Q.1 One-Dimensional Shape Functions

#### Linear (2-node)
```
N₁(ξ) = (1 - ξ) / 2
N₂(ξ) = (1 + ξ) / 2

dN₁/dξ = -1/2
dN₂/dξ = +1/2
```

#### Quadratic (3-node)
```
N₁(ξ) = ξ(ξ - 1) / 2
N₂(ξ) = 1 - ξ²
N₃(ξ) = ξ(ξ + 1) / 2

dN₁/dξ = ξ - 1/2
dN₂/dξ = -2ξ
dN₃/dξ = ξ + 1/2
```

#### Cubic (4-node)
```
N₁(ξ) = 9(1 - ξ)(1/9 - ξ²) / 16
N₂(ξ) = 27(1 - ξ²)(1/3 - ξ) / 16
N₃(ξ) = 27(1 - ξ²)(1/3 + ξ) / 16
N₄(ξ) = 9(1 + ξ)(1/9 - ξ²) / 16
```

#### Hermite Cubic (2-node, C1 continuity)
```
H₁(ξ) = (2 - 3ξ + ξ³) / 4        (displacement at node 1)
H₂(ξ) = L(1 - ξ - ξ² + ξ³) / 8  (slope at node 1)
H₃(ξ) = (2 + 3ξ - ξ³) / 4        (displacement at node 2)
H₄(ξ) = L(-1 - ξ + ξ² + ξ³) / 8 (slope at node 2)

where L = element length, ξ = 2x/L - 1
```

### Q.2 Triangular Shape Functions

#### Linear Triangle (T3)
```
N₁(L₁, L₂, L₃) = L₁
N₂(L₁, L₂, L₃) = L₂
N₃(L₁, L₂, L₃) = L₃

where L₁ + L₂ + L₃ = 1 (area coordinates)

Conversion to Cartesian:
x = L₁x₁ + L₂x₂ + L₃x₃
y = L₁y₁ + L₂y₂ + L₃y₃

Jacobian:
J = | x₂-x₁  y₂-y₁ |
    | x₃-x₁  y₃-y₁ |

Area = det(J)/2 = |J|/2

Derivatives:
∂N/∂x = (1/2A) × [y₂-y₃, y₃-y₁, y₁-y₂]ᵀ
∂N/∂y = (1/2A) × [x₃-x₂, x₁-x₃, x₂-x₁]ᵀ
```

#### Quadratic Triangle (T6)
```
Corner nodes (vertices):
N₁ = L₁(2L₁ - 1)
N₂ = L₂(2L₂ - 1)
N₃ = L₃(2L₃ - 1)

Midside nodes (edge midpoints):
N₄ = 4L₁L₂
N₅ = 4L₂L₃
N₆ = 4L₃L₁

Derivatives:
∂N₁/∂L₁ = 4L₁ - 1
∂N₄/∂L₁ = 4L₂,  ∂N₄/∂L₂ = 4L₁
... etc.
```

### Q.3 Quadrilateral Shape Functions

#### Bilinear Quadrilateral (Q4)
```
N₁(ξ, η) = (1 - ξ)(1 - η) / 4
N₂(ξ, η) = (1 + ξ)(1 - η) / 4
N₃(ξ, η) = (1 + ξ)(1 + η) / 4
N₄(ξ, η) = (1 - ξ)(1 + η) / 4

∂N₁/∂ξ = -(1 - η)/4,  ∂N₁/∂η = -(1 - ξ)/4
∂N₂/∂ξ = +(1 - η)/4,  ∂N₂/∂η = -(1 + ξ)/4
∂N₃/∂ξ = +(1 + η)/4,  ∂N₃/∂η = +(1 + ξ)/4
∂N₄/∂ξ = -(1 + η)/4,  ∂N₄/∂η = +(1 - ξ)/4

Node ordering (counterclockwise):
4----3
|    |
|    |
1----2
```

#### Serendipity Quadrilateral (Q8)
```
Corner nodes (ξᵢ, ηᵢ = ±1):
Nᵢ = (1 + ξξᵢ)(1 + ηηᵢ)(ξξᵢ + ηηᵢ - 1) / 4

Midside nodes on ξ edges (ηᵢ = ±1, ξᵢ = 0):
Nᵢ = (1 - ξ²)(1 + ηηᵢ) / 2

Midside nodes on η edges (ξᵢ = ±1, ηᵢ = 0):
Nᵢ = (1 + ξξᵢ)(1 - η²) / 2

Node ordering:
4--7--3
|     |
8     6
|     |
1--5--2
```

#### Lagrangian Quadrilateral (Q9)
```
Corner nodes (ξᵢ, ηᵢ = ±1):
Nᵢ = ξ(ξ + ξᵢ)η(η + ηᵢ) / 4

Midside nodes on edges:
N₅ = (1 - ξ²)η(η - 1) / 2
N₆ = ξ(ξ + 1)(1 - η²) / 2
N₇ = (1 - ξ²)η(η + 1) / 2
N₈ = ξ(ξ - 1)(1 - η²) / 2

Center node:
N₉ = (1 - ξ²)(1 - η²)

Node ordering:
4--7--3
|     |
8  9  6
|     |
1--5--2
```

### Q.4 Tetrahedral Shape Functions

#### Linear Tetrahedron (T4)
```
N₁ = L₁
N₂ = L₂
N₃ = L₃
N₄ = L₄

where L₁ + L₂ + L₃ + L₄ = 1 (volume coordinates)

Volume = det([x₂-x₁, x₃-x₁, x₄-x₁;
               y₂-y₁, y₃-y₁, y₄-y₁;
               z₂-z₁, z₃-z₁, z₄-z₁]) / 6
```

#### Quadratic Tetrahedron (T10)
```
Corner nodes:
N₁ = L₁(2L₁ - 1)
N₂ = L₂(2L₂ - 1)
N₃ = L₃(2L₃ - 1)
N₄ = L₄(2L₄ - 1)

Edge midpoint nodes:
N₅ = 4L₁L₂   (edge 1-2)
N₆ = 4L₂L₃   (edge 2-3)
N₇ = 4L₃L₁   (edge 3-1)
N₈ = 4L₁L₄   (edge 1-4)
N₉ = 4L₂L₄   (edge 2-4)
N₁₀ = 4L₃L₄  (edge 3-4)
```

### Q.5 Hexahedral Shape Functions

#### Trilinear Hexahedron (H8)
```
Nᵢ(ξ, η, ζ) = (1 + ξξᵢ)(1 + ηηᵢ)(1 + ζζᵢ) / 8

where (ξᵢ, ηᵢ, ζᵢ) = (±1, ±1, ±1) for i = 1, ..., 8

Node ordering (right-hand rule):
Bottom face (ζ = -1):    Top face (ζ = +1):
4----3                   8----7
|    |                   |    |
|    |                   |    |
1----2                   5----6
```

#### Serendipity Hexahedron (H20)
```
Corner nodes (all coordinates ±1):
Nᵢ = (1 + ξξᵢ)(1 + ηηᵢ)(1 + ζζᵢ)(ξξᵢ + ηηᵢ + ζζᵢ - 2) / 8

Edge midpoint nodes (one coordinate = 0, others ±1):
For ξ = 0 edges: Nᵢ = (1 - ξ²)(1 + ηηᵢ)(1 + ζζᵢ) / 4
For η = 0 edges: Nᵢ = (1 + ξξᵢ)(1 - η²)(1 + ζζᵢ) / 4
For ζ = 0 edges: Nᵢ = (1 + ξξᵢ)(1 + ηηᵢ)(1 - ζ²) / 4
```

#### Lagrangian Hexahedron (H27)
```
Uses tensor product of 1D quadratic Lagrange polynomials:

ℓ₁(ξ) = ξ(ξ - 1)/2  (node at ξ = -1)
ℓ₂(ξ) = 1 - ξ²       (node at ξ = 0)
ℓ₃(ξ) = ξ(ξ + 1)/2  (node at ξ = +1)

Shape function for node (i, j, k):
N_{ijk}(ξ, η, ζ) = ℓᵢ(ξ) × ℓⱼ(η) × ℓₖ(ζ)

27 nodes: 8 corners + 12 edge midpoints + 6 face centers + 1 center
```

### Q.6 Implementation

```python
class ShapeFunctions:
    """
    Library of shape functions for finite elements.
    """
    
    @staticmethod
    def t3(L: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Linear triangle shape functions.
        
        Parameters
        ----------
        L : ndarray, shape (3,) or (nip, 3)
            Area coordinates [L1, L2, L3]
        
        Returns
        -------
        N : ndarray, shape (3,) or (nip, 3)
            Shape function values
        dN : ndarray, shape (3, 2) or (nip, 3, 2)
            Shape function derivatives w.r.t. L1, L2
        """
        L = np.atleast_2d(L)
        nip = L.shape[0]
        
        N = L  # N = [L1, L2, L3]
        
        # dN/dL is constant for linear triangle
        dN = np.array([[1, 0],
                       [0, 1],
                       [-1, -1]])
        dN = np.tile(dN, (nip, 1, 1))
        
        return N, dN
    
    @staticmethod
    def t6(L: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Quadratic triangle shape functions.
        """
        L = np.atleast_2d(L)
        nip = L.shape[0]
        L1, L2, L3 = L[:, 0], L[:, 1], L[:, 2]
        
        N = np.zeros((nip, 6))
        N[:, 0] = L1 * (2*L1 - 1)
        N[:, 1] = L2 * (2*L2 - 1)
        N[:, 2] = L3 * (2*L3 - 1)
        N[:, 3] = 4 * L1 * L2
        N[:, 4] = 4 * L2 * L3
        N[:, 5] = 4 * L3 * L1
        
        # Derivatives w.r.t. L1, L2
        dN = np.zeros((nip, 6, 2))
        dN[:, 0, 0] = 4*L1 - 1
        dN[:, 1, 1] = 4*L2 - 1
        dN[:, 2, :] = np.column_stack([-4*L3 + 1, -4*L3 + 1])
        dN[:, 3, 0] = 4*L2
        dN[:, 3, 1] = 4*L1
        dN[:, 4, 0] = -4*L3
        dN[:, 4, 1] = 4*(L3 - L2)
        dN[:, 5, 0] = 4*(L3 - L1)
        dN[:, 5, 1] = -4*L3
        
        return N, dN
    
    @staticmethod
    def q4(xi: float, eta: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Bilinear quadrilateral shape functions.
        """
        N = np.array([
            (1 - xi) * (1 - eta) / 4,
            (1 + xi) * (1 - eta) / 4,
            (1 + xi) * (1 + eta) / 4,
            (1 - xi) * (1 + eta) / 4
        ])
        
        dN = np.array([
            [-(1 - eta)/4, -(1 - xi)/4],
            [+(1 - eta)/4, -(1 + xi)/4],
            [+(1 + eta)/4, +(1 + xi)/4],
            [-(1 + eta)/4, +(1 - xi)/4]
        ])
        
        return N, dN
    
    @staticmethod
    def q8(xi: float, eta: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Serendipity 8-node quadrilateral shape functions.
        """
        N = np.zeros(8)
        dN = np.zeros((8, 2))
        
        # Corner nodes
        xi_i = np.array([-1, 1, 1, -1])
        eta_i = np.array([-1, -1, 1, 1])
        
        for i in range(4):
            xxi = xi * xi_i[i]
            eeta = eta * eta_i[i]
            N[i] = (1 + xxi) * (1 + eeta) * (xxi + eeta - 1) / 4
            dN[i, 0] = xi_i[i] * (1 + eeta) * (2*xxi + eeta) / 4
            dN[i, 1] = eta_i[i] * (1 + xxi) * (xxi + 2*eeta) / 4
        
        # Midside nodes
        N[4] = (1 - xi**2) * (1 - eta) / 2
        N[5] = (1 + xi) * (1 - eta**2) / 2
        N[6] = (1 - xi**2) * (1 + eta) / 2
        N[7] = (1 - xi) * (1 - eta**2) / 2
        
        dN[4, :] = [-xi * (1 - eta), -(1 - xi**2) / 2]
        dN[5, :] = [(1 - eta**2) / 2, -(1 + xi) * eta]
        dN[6, :] = [-xi * (1 + eta), (1 - xi**2) / 2]
        dN[7, :] = [-(1 - eta**2) / 2, -(1 - xi) * eta]
        
        return N, dN
    
    @staticmethod
    def h8(xi: float, eta: float, zeta: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Trilinear hexahedron shape functions.
        """
        xi_i = np.array([-1, 1, 1, -1, -1, 1, 1, -1])
        eta_i = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
        zeta_i = np.array([-1, -1, -1, -1, 1, 1, 1, 1])
        
        N = np.zeros(8)
        dN = np.zeros((8, 3))
        
        for i in range(8):
            N[i] = (1 + xi*xi_i[i]) * (1 + eta*eta_i[i]) * (1 + zeta*zeta_i[i]) / 8
            dN[i, 0] = xi_i[i] * (1 + eta*eta_i[i]) * (1 + zeta*zeta_i[i]) / 8
            dN[i, 1] = eta_i[i] * (1 + xi*xi_i[i]) * (1 + zeta*zeta_i[i]) / 8
            dN[i, 2] = zeta_i[i] * (1 + xi*xi_i[i]) * (1 + eta*eta_i[i]) / 8
        
        return N, dN
```

---

## Appendix R: Strain-Displacement Matrices

### R.1 Two-Dimensional Plane Elements

#### Plane Stress/Strain B-Matrix

For a 2D element with nodal DOFs [u₁, v₁, u₂, v₂, ...], the strain vector is:

```
ε = { εₓₓ }   { ∂u/∂x     }
    { εᵧᵧ } = { ∂v/∂y     }
    { γₓᵧ }   { ∂u/∂y + ∂v/∂x }
```

B-matrix structure:
```
B = [ ∂N₁/∂x    0       ∂N₂/∂x    0      ... ]
    [   0     ∂N₁/∂y     0      ∂N₂/∂y  ... ]
    [ ∂N₁/∂y  ∂N₁/∂x   ∂N₂/∂y  ∂N₂/∂x  ... ]
```

Dimensions: B is (3 × 2n) for n nodes

#### Jacobian and Coordinate Transformation

```
J = [ ∂x/∂ξ  ∂y/∂ξ ] = [ Σ ∂Nᵢ/∂ξ × xᵢ   Σ ∂Nᵢ/∂ξ × yᵢ ]
    [ ∂x/∂η  ∂y/∂η ]   [ Σ ∂Nᵢ/∂η × xᵢ   Σ ∂Nᵢ/∂η × yᵢ ]

J⁻¹ = [ ∂ξ/∂x  ∂η/∂x ]
      [ ∂ξ/∂y  ∂η/∂y ]

∂Nᵢ/∂x = J⁻¹₁₁ × ∂Nᵢ/∂ξ + J⁻¹₁₂ × ∂Nᵢ/∂η
∂Nᵢ/∂y = J⁻¹₂₁ × ∂Nᵢ/∂ξ + J⁻¹₂₂ × ∂Nᵢ/∂η
```

### R.2 Three-Dimensional Solid Elements

#### 3D B-Matrix

For a 3D element with DOFs [u₁, v₁, w₁, u₂, v₂, w₂, ...]:

```
ε = { εₓₓ }   { ∂u/∂x               }
    { εᵧᵧ }   { ∂v/∂y               }
    { εᵧᵧ } = { ∂w/∂z               }
    { γₓᵧ }   { ∂u/∂y + ∂v/∂x       }
    { γᵧᵧ }   { ∂v/∂z + ∂w/∂y       }
    { γₓᵧ }   { ∂u/∂z + ∂w/∂x       }

B = [ ∂N/∂x    0       0     ]
    [   0    ∂N/∂y     0     ]
    [   0      0     ∂N/∂z   ]
    [ ∂N/∂y  ∂N/∂x     0     ]
    [   0    ∂N/∂z   ∂N/∂y   ]
    [ ∂N/∂z    0     ∂N/∂x   ]
```

Dimensions: B is (6 × 3n) for n nodes

### R.3 Beam Elements

#### Euler-Bernoulli Beam (2D)

DOFs: [u₁, v₁, θ₁, u₂, v₂, θ₂]

Axial strain:
```
εₐ = ∂u/∂x = (1/L) × [-1, 0, 0, 1, 0, 0] × d
```

Curvature:
```
κ = ∂²v/∂x² = (1/L²) × [0, 6ξ, L(3ξ-1), 0, -6ξ, L(3ξ+1)] × d
```

where ξ = 2x/L - 1

B-matrix (strain at fiber y from neutral axis):
```
ε(y) = εₐ - y × κ

B = (1/L) × [-1,      -6yξ/L,      -y(3ξ-1),    1,     6yξ/L,      -y(3ξ+1)]
```

### R.4 Shell Elements

#### Kirchhoff Plate (Bending)

DOFs per node: [w, θₓ, θᵧ]

Curvatures:
```
κ = { κₓₓ  }   { -∂²w/∂x²        }
    { κᵧᵧ  } = { -∂²w/∂y²        }
    { κₓᵧ  }   { -2∂²w/∂x∂y      }
```

For DKT element, rotations are related to displacement:
```
θₓ = ∂w/∂y
θᵧ = -∂w/∂x
```

This constraint is enforced discretely at edge midpoints.

#### Reissner-Mindlin Plate (with Shear)

DOFs per node: [w, θₓ, θᵧ]

Curvatures (same as Kirchhoff):
```
κ = { -∂θᵧ/∂x       }
    { ∂θₓ/∂y        }
    { ∂θₓ/∂x - ∂θᵧ/∂y }
```

Transverse shear strains:
```
γ = { γₓᵧ }   { ∂w/∂x - θᵧ }
    { γᵧᵧ } = { ∂w/∂y + θₓ }
```

MITC uses assumed shear strain to avoid locking.

### R.5 Implementation

```python
def compute_B_matrix_plane(dN_dx: np.ndarray) -> np.ndarray:
    """
    Compute B-matrix for plane stress/strain element.
    
    Parameters
    ----------
    dN_dx : ndarray, shape (nnodes, 2)
        Shape function derivatives [∂N/∂x, ∂N/∂y]
    
    Returns
    -------
    B : ndarray, shape (3, 2*nnodes)
        Strain-displacement matrix
    """
    nnodes = dN_dx.shape[0]
    B = np.zeros((3, 2 * nnodes))
    
    for i in range(nnodes):
        B[0, 2*i] = dN_dx[i, 0]      # εxx = ∂u/∂x
        B[1, 2*i+1] = dN_dx[i, 1]    # εyy = ∂v/∂y
        B[2, 2*i] = dN_dx[i, 1]      # γxy = ∂u/∂y + ∂v/∂x
        B[2, 2*i+1] = dN_dx[i, 0]
    
    return B


def compute_B_matrix_3d(dN_dx: np.ndarray) -> np.ndarray:
    """
    Compute B-matrix for 3D solid element.
    
    Parameters
    ----------
    dN_dx : ndarray, shape (nnodes, 3)
        Shape function derivatives [∂N/∂x, ∂N/∂y, ∂N/∂z]
    
    Returns
    -------
    B : ndarray, shape (6, 3*nnodes)
        Strain-displacement matrix
    """
    nnodes = dN_dx.shape[0]
    B = np.zeros((6, 3 * nnodes))
    
    for i in range(nnodes):
        B[0, 3*i] = dN_dx[i, 0]      # εxx
        B[1, 3*i+1] = dN_dx[i, 1]    # εyy
        B[2, 3*i+2] = dN_dx[i, 2]    # εzz
        B[3, 3*i] = dN_dx[i, 1]      # γxy
        B[3, 3*i+1] = dN_dx[i, 0]
        B[4, 3*i+1] = dN_dx[i, 2]    # γyz
        B[4, 3*i+2] = dN_dx[i, 1]
        B[5, 3*i] = dN_dx[i, 2]      # γxz
        B[5, 3*i+2] = dN_dx[i, 0]
    
    return B


def compute_B_matrix_beam2d(
    L: float,
    xi: float,
    y: float = 0.0
) -> np.ndarray:
    """
    Compute B-matrix for 2D Euler-Bernoulli beam.
    
    Parameters
    ----------
    L : float
        Element length
    xi : float
        Natural coordinate (-1 to 1)
    y : float
        Distance from neutral axis (for fiber strain)
    
    Returns
    -------
    B : ndarray, shape (1, 6)
        Strain-displacement matrix
    """
    B = np.zeros(6)
    
    # Axial strain contribution
    B[0] = -1/L
    B[3] = 1/L
    
    # Bending strain contribution (curvature × y)
    B[1] = -y * 6*xi / L**2
    B[2] = -y * (3*xi - 1) / L
    B[4] = y * 6*xi / L**2
    B[5] = -y * (3*xi + 1) / L
    
    return B.reshape(1, 6)
```

---

## Appendix S: Constitutive Matrices

### S.1 Isotropic Linear Elasticity

#### Plane Stress
```
D = E/(1-ν²) × [ 1   ν    0        ]
               [ ν   1    0        ]
               [ 0   0   (1-ν)/2   ]
```

#### Plane Strain
```
D = E/[(1+ν)(1-2ν)] × [ 1-ν    ν      0          ]
                       [  ν   1-ν     0          ]
                       [  0    0    (1-2ν)/2     ]
```

#### Axisymmetric
```
        | εrr |
ε =     | εθθ |
        | εzz |
        | γrz |

D = E/[(1+ν)(1-2ν)] × [ 1-ν    ν      ν      0          ]
                       [  ν   1-ν     ν      0          ]
                       [  ν    ν    1-ν     0          ]
                       [  0    0      0    (1-2ν)/2     ]
```

#### 3D Elasticity
```
D = E/[(1+ν)(1-2ν)] × [ 1-ν    ν      ν      0      0      0      ]
                       [  ν   1-ν     ν      0      0      0      ]
                       [  ν    ν    1-ν     0      0      0      ]
                       [  0    0      0    (1-2ν)/2 0      0      ]
                       [  0    0      0      0    (1-2ν)/2 0      ]
                       [  0    0      0      0      0    (1-2ν)/2 ]
```

### S.2 Orthotropic Elasticity

#### Plane Stress
```
D = [   E₁/(1-ν₁₂ν₂₁)       ν₂₁E₁/(1-ν₁₂ν₂₁)     0    ]
    [ ν₁₂E₂/(1-ν₁₂ν₂₁)      E₂/(1-ν₁₂ν₂₁)       0    ]
    [        0                    0              G₁₂   ]

where ν₂₁ = ν₁₂ × E₂/E₁ (reciprocal relation)
```

#### 3D Orthotropic (9 constants)
```
Independent parameters: E₁, E₂, E₃, ν₁₂, ν₁₃, ν₂₃, G₁₂, G₁₃, G₂₃

Compliance matrix S = D⁻¹:
S = [ 1/E₁    -ν₂₁/E₂  -ν₃₁/E₃   0       0       0     ]
    [ -ν₁₂/E₁  1/E₂    -ν₃₂/E₃   0       0       0     ]
    [ -ν₁₃/E₁ -ν₂₃/E₂   1/E₃     0       0       0     ]
    [    0       0        0     1/G₂₃    0       0     ]
    [    0       0        0       0     1/G₁₃    0     ]
    [    0       0        0       0       0     1/G₁₂  ]
```

### S.3 Plate/Shell Constitutive Relations

#### Kirchhoff Plate Bending
```
M = { Mxx }   D_b = Et³/[12(1-ν²)] × [ 1   ν   0       ]
    { Myy } =                         [ ν   1   0       ] × κ
    { Mxy }                           [ 0   0  (1-ν)/2  ]
```

#### Reissner-Mindlin with Shear
```
Bending:  M = D_b × κ    (same as above)
Shear:    Q = D_s × γ

D_s = κ × G × t × [ 1  0 ]
                   [ 0  1 ]

where κ = 5/6 (shear correction factor)
```

### S.4 Implementation

```python
def D_plane_stress(E: float, nu: float) -> np.ndarray:
    """
    Constitutive matrix for plane stress.
    
    Parameters
    ----------
    E : float
        Young's modulus
    nu : float
        Poisson's ratio
    
    Returns
    -------
    D : ndarray, shape (3, 3)
        Constitutive matrix
    """
    factor = E / (1 - nu**2)
    D = factor * np.array([
        [1,   nu,       0           ],
        [nu,  1,        0           ],
        [0,   0,        (1-nu)/2    ]
    ])
    return D


def D_plane_strain(E: float, nu: float) -> np.ndarray:
    """
    Constitutive matrix for plane strain.
    """
    factor = E / ((1 + nu) * (1 - 2*nu))
    D = factor * np.array([
        [1-nu,  nu,       0           ],
        [nu,    1-nu,     0           ],
        [0,     0,        (1-2*nu)/2  ]
    ])
    return D


def D_3d_isotropic(E: float, nu: float) -> np.ndarray:
    """
    Constitutive matrix for 3D isotropic elasticity.
    """
    factor = E / ((1 + nu) * (1 - 2*nu))
    c1 = 1 - nu
    c2 = nu
    c3 = (1 - 2*nu) / 2
    
    D = factor * np.array([
        [c1, c2, c2,  0,  0,  0],
        [c2, c1, c2,  0,  0,  0],
        [c2, c2, c1,  0,  0,  0],
        [ 0,  0,  0, c3,  0,  0],
        [ 0,  0,  0,  0, c3,  0],
        [ 0,  0,  0,  0,  0, c3]
    ])
    return D


def D_plate_bending(E: float, nu: float, t: float) -> np.ndarray:
    """
    Constitutive matrix for Kirchhoff plate bending.
    
    Parameters
    ----------
    E : float
        Young's modulus
    nu : float
        Poisson's ratio
    t : float
        Plate thickness
    
    Returns
    -------
    D : ndarray, shape (3, 3)
        Moment-curvature constitutive matrix
    """
    factor = E * t**3 / (12 * (1 - nu**2))
    D = factor * np.array([
        [1,   nu,       0           ],
        [nu,  1,        0           ],
        [0,   0,        (1-nu)/2    ]
    ])
    return D


def D_plate_shear(G: float, t: float, kappa: float = 5/6) -> np.ndarray:
    """
    Constitutive matrix for Reissner-Mindlin shear.
    
    Parameters
    ----------
    G : float
        Shear modulus
    t : float
        Plate thickness
    kappa : float
        Shear correction factor (default 5/6)
    
    Returns
    -------
    D : ndarray, shape (2, 2)
        Shear constitutive matrix
    """
    return kappa * G * t * np.eye(2)
```

---

## Document Statistics

| Metric | Value |
|--------|-------|
| Total lines | ~10,500 |
| Total tasks | 247 |
| Element tasks | 89 |
| Material tasks | 78 |
| Analysis tasks | 52 |
| Testing tasks | 28 |
| Estimated effort | 18-24 months |
| Priority CRITICAL | 43 tasks |
| Priority HIGH | 87 tasks |
| Priority MEDIUM | 72 tasks |
| Priority LOW | 45 tasks |

---

**END OF DOCUMENT**

*Document generated: 2026-04-05*
*Author: GitHub Copilot CLI*
*Repository: github.com/adzetto/femlabpy*
