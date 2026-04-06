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

This chapter develops the mathematical foundations for plasticity models in `femlabpy.materials`. We derive the von Mises and Drucker-Prager yield criteria, establish the return mapping algorithms, and connect them to the implementation.

## 1. Stress Invariants and Decomposition

### 1.1 Stress tensor decomposition

Any stress tensor $\boldsymbol{\sigma}$ can be decomposed into volumetric (hydrostatic) and deviatoric parts:

$$
\boldsymbol{\sigma} = \underbrace{\sigma_m \mathbf{I}}_{\text{volumetric}} + \underbrace{\mathbf{s}}_{\text{deviatoric}}
$$

where the mean (hydrostatic) stress is:

$$
\sigma_m = \frac{1}{3} \text{tr}(\boldsymbol{\sigma}) = \frac{1}{3}(\sigma_{xx} + \sigma_{yy} + \sigma_{zz})
$$

and the deviatoric stress tensor is:

$$
\mathbf{s} = \boldsymbol{\sigma} - \sigma_m \mathbf{I}
$$

In Voigt notation for 3D: $\mathbf{S} = [\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{zx}]^T$, the deviatoric components are:

$$
s_{ii} = \sigma_{ii} - \sigma_m, \quad s_{ij} = \sigma_{ij} \text{ for } i \neq j
$$

### 1.2 Principal stress invariants

The three invariants of the stress tensor are:

$$
I_1 = \text{tr}(\boldsymbol{\sigma}) = \sigma_{xx} + \sigma_{yy} + \sigma_{zz}
$$

$$
I_2 = \frac{1}{2}\left[ (\text{tr}\boldsymbol{\sigma})^2 - \text{tr}(\boldsymbol{\sigma}^2) \right]
$$

$$
I_3 = \det(\boldsymbol{\sigma})
$$

The deviatoric invariants $J_2$ and $J_3$ are defined as:

$$
J_2 = \frac{1}{2} \mathbf{s} : \mathbf{s} = \frac{1}{2} s_{ij} s_{ij}
$$

$$
J_3 = \det(\mathbf{s})
$$

Expanding $J_2$ in terms of Voigt components:

$$
J_2 = \frac{1}{2}\left( s_{xx}^2 + s_{yy}^2 + s_{zz}^2 \right) + s_{xy}^2 + s_{yz}^2 + s_{zx}^2
$$

### 1.3 Von Mises equivalent stress

The von Mises equivalent stress (or effective stress) is defined as:

$$
\sigma_{eq} = \sqrt{3 J_2} = \sqrt{\frac{3}{2} \mathbf{s} : \mathbf{s}}
$$

This can be written explicitly as:

$$
\sigma_{eq} = \sqrt{\frac{1}{2}\left[ (\sigma_{xx}-\sigma_{yy})^2 + (\sigma_{yy}-\sigma_{zz})^2 + (\sigma_{zz}-\sigma_{xx})^2 \right] + 3(\sigma_{xy}^2 + \sigma_{yz}^2 + \sigma_{zx}^2)}
$$

For plane stress ($\sigma_{zz} = \sigma_{yz} = \sigma_{zx} = 0$):

$$
\sigma_{eq} = \sqrt{\sigma_{xx}^2 - \sigma_{xx}\sigma_{yy} + \sigma_{yy}^2 + 3\sigma_{xy}^2}
$$

This is implemented in `eqstress(S)`.

---

## 2. Von Mises Yield Criterion

### 2.1 Yield function

The von Mises yield criterion states that yielding occurs when the equivalent stress equals the yield stress:

$$
f(\boldsymbol{\sigma}) = \sigma_{eq} - \sigma_y = \sqrt{3 J_2} - \sigma_y = 0
$$

This defines a cylindrical yield surface in principal stress space, with axis along the hydrostatic line $\sigma_1 = \sigma_2 = \sigma_3$.

### 2.2 Isotropic hardening

With linear isotropic hardening, the yield stress evolves with accumulated plastic strain:

$$
\sigma_y = \sigma_{y0} + H \bar{\varepsilon}^p
$$

where:
- $\sigma_{y0}$ is the initial yield stress
- $H$ is the hardening modulus
- $\bar{\varepsilon}^p$ is the equivalent plastic strain

The increment of equivalent plastic strain is related to the plastic multiplier $\Delta\lambda$ by:

$$
\Delta \bar{\varepsilon}^p = \Delta\lambda
$$

### 2.3 Associative flow rule

The plastic strain rate follows the normality rule (associated flow):

$$
\dot{\boldsymbol{\varepsilon}}^p = \dot{\lambda} \frac{\partial f}{\partial \boldsymbol{\sigma}} = \dot{\lambda} \frac{\partial \sigma_{eq}}{\partial \boldsymbol{\sigma}} = \dot{\lambda} \frac{3 \mathbf{s}}{2 \sigma_{eq}}
$$

The direction tensor (normal to the yield surface) is:

$$
\mathbf{n} = \frac{\partial f}{\partial \boldsymbol{\sigma}} = \frac{3 \mathbf{s}}{2 \sigma_{eq}}
$$

### 2.4 Elastic-plastic tangent modulus

The consistent (algorithmic) tangent for 3D is:

$$
\mathbf{D}^{ep} = \mathbf{D} - \frac{(\mathbf{D} : \mathbf{n}) \otimes (\mathbf{n} : \mathbf{D})}{\mathbf{n} : \mathbf{D} : \mathbf{n} + H}
$$

For radial return in 3D von Mises plasticity, this simplifies because the elastic stiffness $\mathbf{D}$ can be split into deviatoric and volumetric parts.

---

## 3. Return Mapping Algorithm

### 3.1 Elastic predictor

Given the strain increment $\Delta\boldsymbol{\varepsilon}$, compute the trial (elastic) stress:

$$
\boldsymbol{\sigma}^{trial} = \boldsymbol{\sigma}_n + \mathbf{D} : \Delta\boldsymbol{\varepsilon}
$$

where $\boldsymbol{\sigma}_n$ is the stress at the end of the previous increment.

### 3.2 Yield check

Evaluate the yield function at the trial state:

$$
f^{trial} = \sigma_{eq}^{trial} - \sigma_{y,n}
$$

If $f^{trial} \leq 0$, the step is elastic: $\boldsymbol{\sigma}_{n+1} = \boldsymbol{\sigma}^{trial}$.

If $f^{trial} > 0$, plastic correction is needed.

### 3.3 Plastic corrector (3D radial return)

For J2 plasticity in 3D, the return mapping has a closed-form solution. The deviatoric trial stress is:

$$
\mathbf{s}^{trial} = \boldsymbol{\sigma}^{trial} - \sigma_m^{trial} \mathbf{I}
$$

The plastic multiplier is found by enforcing the consistency condition:

$$
f(\boldsymbol{\sigma}_{n+1}) = 0 \implies \sigma_{eq,n+1} = \sigma_{y,n+1}
$$

This gives:

$$
\Delta\lambda = \frac{f^{trial}}{3G + H}
$$

where $G = \frac{E}{2(1+\nu)}$ is the shear modulus.

The corrected deviatoric stress is:

$$
\mathbf{s}_{n+1} = \left(1 - \frac{3G \Delta\lambda}{\sigma_{eq}^{trial}}\right) \mathbf{s}^{trial}
$$

The updated stress is:

$$
\boldsymbol{\sigma}_{n+1} = \mathbf{s}_{n+1} + \sigma_m^{trial} \mathbf{I}
$$

### 3.4 Plane stress return mapping

Under plane stress constraints ($\sigma_{zz} = 0$), the radial return is no longer closed-form. The constraint must be enforced iteratively.

The residual system involves finding $\Delta\lambda$ such that:

$$
\begin{cases}
f = \sigma_{eq} - \sigma_y = 0 \\
\sigma_{zz} = 0
\end{cases}
$$

The implementation uses Newton-Raphson iteration on the scalar consistency equation. Functions `yieldvm(S, G, dL, Sy)` and `dyieldvm(S, G, dL, Sy)` evaluate the residual and its derivative.

The iteration proceeds:

$$
\Delta\lambda^{(k+1)} = \Delta\lambda^{(k)} - \frac{f(\Delta\lambda^{(k)})}{\frac{\partial f}{\partial \Delta\lambda}\big|_{\Delta\lambda^{(k)}}}
$$

until $|f| < \text{tol}$.

---

## 4. Drucker-Prager Plasticity

The Drucker-Prager criterion models pressure-dependent yielding, relevant for soils, rocks, and concrete.

### 4.1 Yield function

$$
f(\boldsymbol{\sigma}) = \sigma_{eq} + \phi \sigma_m - \sigma_y = \sqrt{3J_2} + \phi \cdot \frac{I_1}{3} - \sigma_y = 0
$$

where $\phi$ is the friction parameter controlling pressure sensitivity.

When $\phi = 0$, this reduces to von Mises plasticity.

### 4.2 Geometric interpretation

In principal stress space, the Drucker-Prager surface is a cone with:
- Apex on the hydrostatic axis
- Cone angle determined by $\phi$
- Von Mises cylinder as the limiting case ($\phi \to 0$)

### 4.3 Flow rule and derivatives

The gradient of the yield function:

$$
\frac{\partial f}{\partial \boldsymbol{\sigma}} = \frac{3 \mathbf{s}}{2 \sigma_{eq}} + \frac{\phi}{3} \mathbf{I}
$$

This non-deviatoric term means plastic dilation (volume change) accompanies yielding.

### 4.4 Drucker-Prager return mapping

The `stressdp` function implements a coupled Newton iteration solving for both stress correction $\Delta\mathbf{S}$ and plastic multiplier $\Delta\lambda$ simultaneously.

The residual vector is:

$$
\mathbf{R} = \Delta\boldsymbol{\varepsilon} - \mathbf{C} : (\Delta\boldsymbol{\sigma} + \delta\boldsymbol{\sigma}) - \Delta\lambda \frac{\partial f}{\partial \boldsymbol{\sigma}}
$$

The $4 \times 4$ tangent system (for 3 stress components + multiplier) is solved at each iteration:

$$
\begin{bmatrix} \mathbf{C} + \Delta\lambda \frac{\partial^2 f}{\partial \boldsymbol{\sigma}^2} & \frac{\partial f}{\partial \boldsymbol{\sigma}} \\ \left(\frac{\partial f}{\partial \boldsymbol{\sigma}}\right)^T & -H \end{bmatrix} \begin{bmatrix} \delta(\Delta\boldsymbol{\sigma}) \\ \delta(\Delta\lambda) \end{bmatrix} = -\begin{bmatrix} \mathbf{R} \\ f \end{bmatrix}
$$

---

## 5. Material Row Layout

The plasticity routines expect material parameters in a flat array:

| Index | Parameter | Description |
|-------|-----------|-------------|
| 0     | $E$       | Young's modulus |
| 1     | $\nu$     | Poisson's ratio |
| 2     | $\sigma_{y0}$ | Initial yield stress |
| 3     | $H$       | Hardening modulus |
| 4     | $\phi$    | Friction parameter (Drucker-Prager only) |

---

## 6. Implementation Summary

### Module structure

| Function | Description |
|----------|-------------|
| `devstress(S)` | Splits stress into deviatoric and mean parts |
| `eqstress(S)` | Computes von Mises equivalent stress |
| `yieldvm(S, G, dL, Sy)` | Evaluates plane-stress consistency residual |
| `dyieldvm(S, G, dL, Sy)` | Derivative of consistency residual |
| `stressvm(S, G, Sy)` | Plane-stress von Mises return mapping |
| `stressdp(S, G, Sy0, dE, dS)` | Drucker-Prager return mapping |

### Design principles

1. **Pure functions**: Invariant helpers do not mutate inputs
2. **Explicit state**: Trial stress, plastic multiplier, and corrected stress are separate variables
3. **Local solvers**: Each material call is self-contained; no global state caching
4. **Modular tangent**: Consistent tangent computed at solver level, not material level

### Source file locations

```
femlabpy/
├── materials/
│   ├── __init__.py      # Public exports, devstres alias
│   ├── invariants.py    # devstress, eqstress
│   └── plasticity.py    # yieldvm, dyieldvm, stressvm, stressdp
```
