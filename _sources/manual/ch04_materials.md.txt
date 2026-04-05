# Chapter 4: Material Models

`femlabpy` keeps material behavior intentionally compact at the API level. Most routines work with a small material row `G`, a stress vector `S`, and, when plasticity is active, a few scalar history variables. That makes the data flow easy to trace from element routines to Gauss-point updates and then back into the global solver.

This chapter focuses on three things:

1. How the library represents elastic and plastic material data.
2. What the public material helper functions do.
3. How the local return-mapping variables fit together during an update.

## 4.1 Material Data Layout

The material models in this repository are table driven. Instead of passing a large object tree, the element routines expect a compact NumPy row. That choice matches the legacy FemLab style and keeps the kernels easy to vectorize.

For linear elasticity, the important entries are:

- `E`: Young's modulus.
- `nu`: Poisson's ratio.

For the plasticity helpers, the same row is extended with:

- `Sy0`: initial yield stress.
- `H`: isotropic hardening modulus.
- `phi`: Drucker-Prager friction parameter, where that model is used.

The stress vectors follow the same compact convention:

- Plane stress and plane strain routines typically use `[sxx, syy, sxy]`.
- Three-dimensional invariant helpers accept `[sxx, syy, szz, syz, sxz, sxy]`.

That layout matters because most helper functions in `src/femlabpy/materials` are thin transformations on these arrays. They do not hide the data; they expose it in a form that the element and solver code can reuse directly.

### Public helper functions

The `femlabpy.materials` package re-exports a small set of helpers:

- `devstress`
- `devstres` as a compatibility alias for `devstress`
- `eqstress`
- `yieldvm`
- `dyieldvm`
- `stressvm`
- `stressdp`

The first two are invariant utilities. The next three implement the von Mises consistency equation and its return mapping. The last one handles Drucker-Prager plastic correction.

## 4.2 Linear Elasticity

For isotropic linear elasticity, the constitutive law is written in matrix form as:

$$ \sigma = \mathbf{D} \epsilon $$

The exact form of `D` depends on the kinematic assumption.

### Plane Stress

Plane stress is used when the out-of-plane stress components are neglected:

$$ \sigma_z = \tau_{xz} = \tau_{yz} = 0 $$

The corresponding material matrix is:

$$ \mathbf{D} = \frac{E}{1 - \nu^2}
\begin{bmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \frac{1 - \nu}{2}
\end{bmatrix} $$

This is the standard choice for thin plates and membrane-like 2D idealizations. In the codebase, it is the form most closely aligned with the plane-stress plasticity update.

### Plane Strain

Plane strain assumes the out-of-plane strains vanish:

$$ \epsilon_z = \gamma_{xz} = \gamma_{yz} = 0 $$

The constitutive matrix becomes:

$$ \mathbf{D} = \frac{E}{(1 + \nu)(1 - 2\nu)}
\begin{bmatrix}
1-\nu & \nu & 0 \\
\nu & 1-\nu & 0 \\
0 & 0 & \frac{1 - 2\nu}{2}
\end{bmatrix} $$

This is the form used for thick sections, long bodies, or any case where the out-of-plane deformation is strongly constrained.

### Why the chapter still matters for plasticity

Even when a material is not purely elastic, the elastic matrix is still the backbone of the update. It supplies the trial stress, the compliance matrix, and the local stiffness that the solver uses when it linearizes the constitutive response.

## 4.3 Stress Invariants

Plasticity models in this repository are built around invariants rather than raw component-by-component checks. That is the right choice for isotropic yield criteria because the yield surface depends on the stress state as a whole, not on one stress component in isolation.

### Deviatoric stress

The deviatoric helper follows the compact legacy convention used by the codebase: it strips out the mean-stress contribution in the form expected by the plane and 3D helper routines. In the code, `devstress(S)` returns two values:

- the deviatoric stress vector
- the mean stress

That split is useful because von Mises plasticity only cares about the deviatoric part, while Drucker-Prager also depends on the mean stress.

### Equivalent stress

`eqstress(S)` computes the von Mises equivalent stress from the stress vector. This is the scalar quantity used to compare the current state against the yield stress.

For a 2D stress vector `[sxx, syy, sxy]`, the helper evaluates the standard equivalent-stress expression and returns a single floating-point value. For the 3D form, it uses the six-component tensor convention and includes the shear terms explicitly.

### Why this is separated from the plasticity update

The code keeps invariant calculations in a separate module because they are reused by multiple models:

- `stressvm` uses equivalent-stress logic to drive the von Mises correction.
- `stressdp` uses both deviatoric and mean stress to build the Drucker-Prager residual.

This separation keeps the constitutive code readable and avoids duplicating the same tensor algebra in multiple solvers.

## 4.4 Von Mises Elastoplasticity

The von Mises implementation in `src/femlabpy/materials/plasticity.py` is a legacy plane-stress return-mapping routine. It is built around a single scalar unknown: the plastic multiplier increment `dL`.

### Yield function

`yieldvm(S, G, dL, Sy)` evaluates the return-mapping consistency equation. Conceptually, it answers a yes/no question:

- `f <= 0`: the state is admissible.
- `f > 0`: the trial stress lies outside the yield surface and must be projected back.

The helper is not just a classifier. It is the scalar residual used by Newton's method to solve for `dL`.

### Derivative of the yield function

`dyieldvm(S, G, dL, Sy)` returns the analytical derivative of the consistency equation with respect to `dL`.

That derivative is important because the local solve is not a generic black-box root finder. The return mapping converges faster and more reliably when the derivative is computed directly from the constitutive equations.

### Stress return mapping

`stressvm(S, G, Sy)` is the actual correction step. It performs three jobs:

1. Start from a trial stress state.
2. Solve for `dL` with a scalar Newton loop.
3. Use the converged multiplier to update the stress components.

The code first evaluates the yield residual at `dL = 0`. If the state is already admissible, the loop exits almost immediately. If not, it repeatedly updates `dL` with:

$$ dL \leftarrow dL - \frac{f}{f'} $$

until the residual is small enough.

The main pieces of data carried through the loop are:

- `stress`: the trial stress, stored as a flat vector.
- `material`: the row containing `E`, `nu`, and `H`.
- `Sy`: the current yield stress, which is increased by `H * dL` after convergence.
- `dL`: the accumulated plastic multiplier increment.
- `f`: the current consistency residual.
- `df`: the derivative of that residual.

After the multiplier converges, the routine computes the corrected stress components and returns:

- the updated stress vector as a column array
- the converged scalar `dL`

### What the helper functions mean in practice

This is the part that usually matters when you are tracing a failing element step:

- `yieldvm` tells you whether the trial state is inside or outside the yield surface.
- `dyieldvm` tells you how steep that scalar residual is.
- `stressvm` applies the correction and gives you the projected stress that the element should store at the integration point.

If a plastic step does not converge, the first things to inspect are the trial stress, the hardening modulus, and the sign/size of the derivative returned by `dyieldvm`.

## 4.5 Drucker-Prager Plasticity

`stressdp(S, G, Sy0, dE, dS)` implements a Drucker-Prager correction for pressure-sensitive plasticity. Compared with the von Mises update, this one carries more local state because the yield function depends on both deviatoric and mean stress, and because the local Newton system solves for more than one unknown at a time.

### Inputs used by the local solve

The function works with these pieces of data:

- `S`: the trial stress state.
- `G`: the material row containing `E`, `nu`, `Sy0`, `H`, and `phi`.
- `Sy0`: the yield stress at the start of the increment.
- `dE`: the strain increment at the integration point.
- `dS`: the elastic trial stress increment.

That split is deliberate. The function does not recompute the whole finite element step. It only corrects the local constitutive state once the element kernel has already formed the trial increment.

### Residuals and tangent data

The local solve uses two residuals:

- `f`: the Drucker-Prager yield residual.
- `R`: the constitutive consistency residual for the stress update.

It also builds several intermediate quantities that are worth knowing when debugging:

- `Sd`: deviatoric stress.
- `Sm`: mean stress.
- `Seq`: equivalent stress.
- `df`: gradient of the yield function.
- `d2f`: curvature of the yield function.
- `deltaS`: cumulative stress correction.
- `dL`: plastic multiplier increment.

The code assembles the local Jacobian with `np.block`, which makes the structure explicit:

```python
tangent = np.block([[C + dL * d2f, df], [df.T, np.array([[-H]], dtype=float)]])
delta = np.linalg.solve(tangent, np.vstack([R, [[-f]]]))
```

That block matrix is the local linearization of the coupled stress and multiplier system. It is the same idea as a global Newton step, just restricted to one integration point.

### Why this is a coupled solve

Unlike the scalar von Mises update, Drucker-Prager does not reduce cleanly to a single multiplier equation. Pressure sensitivity makes the correction depend on the current stress direction and the yield-surface geometry. The code therefore iterates on both the stress correction and the multiplier until both residuals are small.

The routine returns:

- the corrected stress vector
- the converged plastic multiplier increment

## 4.6 Worked Example

The following sketch shows how the public helpers fit together when a trial stress exceeds yield:

```python
import numpy as np
from femlabpy.materials import yieldvm, stressvm

# Material row: [E, nu, Sy0, H]
G = np.array([200000.0, 0.3, 250.0, 2000.0])
Sy = G[2]

# Trial stress in plane form: [sxx, syy, sxy]
S_trial = np.array([400.0, 0.0, 0.0])

f_trial = yieldvm(S_trial, G, 0.0, Sy)
if f_trial > 0:
    S_corr, dL = stressvm(S_trial, G, Sy)
```

In a real element loop, the returned stress would feed the internal-force vector, and the updated tangent would be used by the global Newton iteration. The key point is that the material update is local, but its effect is global.

## 4.7 Reading the Source

The implementation is small enough to read directly, and that is the best way to understand the code path:

- [`src/femlabpy/materials/__init__.py`](../../src/femlabpy/materials/__init__.py) shows the public exports.
- [`src/femlabpy/materials/invariants.py`](../../src/femlabpy/materials/invariants.py) contains the invariant helpers.
- [`src/femlabpy/materials/plasticity.py`](../../src/femlabpy/materials/plasticity.py) contains the return-mapping routines.

If you are tracing a solver failure, start from the element kernel, inspect the trial stress, and then follow the values passed into `yieldvm`, `stressvm`, or `stressdp`. That is usually enough to locate whether the issue is in the constitutive law, the integration step, or the global assembly.
