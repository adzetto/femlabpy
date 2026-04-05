---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Constitutive Models and Plasticity

This chapter tracks the current `femlabpy.materials` implementation. The important idea is not just the formulas, but how the code carries state through a local update. The module uses small, explicit helper functions for invariants and return mapping, and each solver call works on copied arrays so the trial state, corrected state, and plastic increment stay separate.

## 1. Elasticity and invariant helpers

The public material API is built from a few small functions:

`devstress`, `devstres`, `eqstress`, `yieldvm`, `dyieldvm`, `stressvm`, and `stressdp`

`devstres` is only an alias for `devstress`. The helpers in `invariants.py` do not mutate their inputs. They reshape the stress vector, compute the deviatoric part, and return the mean stress or equivalent stress as needed.

### 1.1 Stress splitting

`devstress(S)` accepts either plane-stress style input `[sxx, syy, sxy]` or 3D input `[sxx, syy, szz, syz, sxz, sxy]`. The function copies the data, subtracts the hydrostatic mean from the normal components, and returns:

1. The deviatoric stress vector.
2. The mean stress as a scalar.

That split is reused by the Drucker-Prager return mapping.

`eqstress(S)` computes the von Mises equivalent stress for both the 2D and 3D layouts. In the code, it is a scalar post-processing helper, not a stateful object.

### 1.2 Material row layout

The plasticity routines expect a flat material row. The code reads:

`[E, nu, Sy0, H, ...]`

for von Mises plasticity, and

`[E, nu, Sy0, H, phi]`

for Drucker-Prager. The functions only use the entries they need, so the state passed in by the caller must stay consistent with those conventions.

## 2. Plane-stress von Mises plasticity

The plane-stress path is implemented by `yieldvm`, `dyieldvm`, and `stressvm`. These functions are tightly coupled: `yieldvm` defines the scalar consistency equation, `dyieldvm` gives its derivative, and `stressvm` runs the Newton loop that updates the plastic multiplier.

### 2.1 Consistency residual

`yieldvm(S, G, dL, Sy)` evaluates the plane-stress return-mapping residual for a candidate plastic increment `dL`. The function unpacks `E`, `nu`, and `H`, builds the transformed moduli `E1` and `E2`, and then evaluates the scalar residual

`f(dL) = 0`

for the trial stress vector. The routine does not update the stress state itself. It only answers the question: "for this stress and this plastic multiplier, are we on the yield surface yet?"

### 2.2 Analytical derivative

`dyieldvm(S, G, dL, Sy)` differentiates the same residual with respect to `dL`. That is what makes the Newton loop in `stressvm` stable and compact. The derivative is computed directly from the same `E1`, `E2`, and stress terms used by `yieldvm`, so the two functions stay consistent.

### 2.3 Return mapping loop

`stressvm(S, G, Sy)` is where the state actually changes. The implementation:

1. Copies the trial stress into a local `stress` array.
2. Sets `dL = 0.0`.
3. Evaluates the residual with `yieldvm`.
4. Repeats Newton updates until `abs(f) <= 1.0e-6`.
5. Updates the yield stress with `Sy = Sy + H * dL`.
6. Projects the stress components back onto the plane-stress yield surface.

That loop is the whole local state machine. The function returns the corrected stress as a column vector and the converged plastic increment as a float. Nothing is stored globally.

The key design point is that plane-stress plasticity is not handled by a closed-form radial return in the same way as the 3D model. The `sigma_zz = 0` constraint forces the code to iterate on the transformed residual instead.

## 3. Drucker-Prager plasticity

`stressdp(S, G, Sy0, dE, dS)` handles the pressure-dependent model. This function is more explicit about state because it keeps track of the stress correction `deltaS`, the plastic multiplier `dL`, and the local residual vector `R` at the same time.

### 3.1 Local state variables

The routine begins by reshaping the inputs into column vectors:

`stress`, `dE`, and `dS`

It then builds the elastic compliance-like matrix `C` from `E` and `nu`, splits the trial stress into deviatoric and mean parts with `devstress`, and evaluates the equivalent stress with `eqstress`. The yield function is initialized as

`f = Seq + phi * Sm - Sy0`

where `phi` is the friction parameter stored in the material row.

### 3.2 Newton system

The return mapping is driven by a small local Newton solve. At each iteration the code computes:

1. `d2f`, the second derivative of the yield function.
2. A 4 x 4 tangent matrix built from `C`, `df`, `d2f`, and `H`.
3. The Newton correction `delta` from `np.linalg.solve`.

The state update is explicit:

```python
deltaS += delta[0:3]
dL += float(delta[3, 0])
```

Then the code recomputes `Sd`, `Sm`, `Seq`, `Sy`, `f`, `df`, and the residual

`R = dE - C @ (dS + deltaS) - dL * df`

That last line is the main difference from the simpler von Mises path. The Drucker-Prager routine solves for both the stress correction and the multiplier together, so the code keeps them as separate evolving variables until convergence.

### 3.3 What the function returns

`stressdp` returns the corrected stress vector and the converged plastic increment. Like the other material helpers, it does not cache any history internally. The caller must supply the current state every time and must store the returned result if it needs the updated stress for the next increment.

## 4. Consistent tangent in the solver flow

The material module stops at the local stress update. It does not expose a dedicated tangent helper. That is a deliberate boundary: the return map is the state update, while the consistent tangent belongs to the global Newton solver or a higher-level wrapper that differentiates the local update.

In practice, that means `stressvm` and `stressdp` are the functions you call inside a constitutive loop, and the solver decides whether it needs the algorithmic modulus for convergence control. The important point for readers is that the local routines already carry the full history needed for a single increment:

1. Trial state comes in through the function arguments.
2. The local Newton loop updates the plastic multiplier and corrected stress.
3. The corrected state is returned to the caller.
4. Any tangent used by the global solve is derived from that returned state.

## 5. How to read the module

The shortest path through the source is:

1. `invariants.py` for `devstress` and `eqstress`.
2. `plasticity.py` for the plane-stress von Mises and Drucker-Prager updates.
3. `materials/__init__.py` for the public export names and the `devstres` alias.

That split is intentional. The invariant helpers are reusable primitives, while the plasticity routines are local solvers with clear input-output boundaries. If you extend the module, keep that separation: pure invariants first, state updates second.
