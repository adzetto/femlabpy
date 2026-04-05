---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Structural Dynamics and Time Integration

This chapter provides a rigorous mathematical foundation for the solver side of `femlabpy`: mass and damping models, time integration algorithms (Newmark-$\beta$, HHT-$\alpha$, Central Difference), and frequency response analysis.

The semidiscrete equation of motion for a linear structural system is given by the second-order ordinary differential equation:

$$
\mathbf{M}\ddot{\mathbf{u}}(t) + \mathbf{C}\dot{\mathbf{u}}(t) + \mathbf{K}\mathbf{u}(t) = \mathbf{p}(t)
$$

where $\mathbf{M}$, $\mathbf{C}$, and $\mathbf{K}$ are the global mass, damping, and stiffness matrices respectively; $\mathbf{u}(t)$, $\dot{\mathbf{u}}(t)$, and $\ddot{\mathbf{u}}(t)$ are the displacement, velocity, and acceleration vectors; and $\mathbf{p}(t)$ is the external load vector.

## 1. Mass and Damping Models

### 1.1. Mass Matrices (Consistent vs. Lumped)

**Consistent Mass Matrix ($\mathbf{M}_c$):**
Derived using the same shape functions $\mathbf{N}(\mathbf{x})$ used for the stiffness matrix. It is full (non-diagonal) and kinematically coupled:

$$
\mathbf{M}_c = \int_{\Omega^e} \rho \mathbf{N}^T \mathbf{N} \, d\Omega
$$

**Lumped Mass Matrix ($\mathbf{M}_L$):**
A diagonal matrix that decouples inertial forces. It is strictly required for explicit time integration (e.g., Central Difference). Lumping is typically achieved via row-summation or the HRZ (Hinton-Rock-Zienkiewicz) scaling method:

$$
M_{L, ii} = \alpha \int_{\Omega^e} \rho N_i^2 \, d\Omega \quad \text{where} \quad \alpha = \frac{\int \rho \, d\Omega}{\sum \int \rho N_j^2 \, d\Omega}
$$

### 1.2. Rayleigh Damping

Rayleigh damping constructs the global damping matrix as a linear combination of the mass and stiffness matrices:

$$
\mathbf{C} = \alpha \mathbf{M} + \beta \mathbf{K}
$$

The coefficients $\alpha$ (mass-proportional) and $\beta$ (stiffness-proportional) are determined by specifying desired damping ratios $\zeta_1$ and $\zeta_2$ at two target natural frequencies $\omega_1$ and $\omega_2$:

$$
\begin{bmatrix} \alpha \\ \beta \end{bmatrix} = \frac{2}{\omega_2^2 - \omega_1^2} \begin{bmatrix} \omega_2 & -\omega_1 \\ -1/\omega_2 & 1/\omega_1 \end{bmatrix} \begin{bmatrix} \omega_1 \zeta_1 \\ \omega_2 \zeta_2 \end{bmatrix}
$$

## 2. Implicit Time Integration: Newmark-$\beta$ Method

The Newmark method relies on the following finite difference expansions for displacement and velocity:

$$
\mathbf{u}_{n+1} = \mathbf{u}_n + \Delta t \dot{\mathbf{u}}_n + \frac{\Delta t^2}{2} \left[ (1 - 2\beta)\ddot{\mathbf{u}}_n + 2\beta \ddot{\mathbf{u}}_{n+1} \right]
$$

$$
\dot{\mathbf{u}}_{n+1} = \dot{\mathbf{u}}_n + \Delta t \left[ (1 - \gamma)\ddot{\mathbf{u}}_n + \gamma \ddot{\mathbf{u}}_{n+1} \right]
$$

By substituting these into the equation of motion at time $t_{n+1}$, we form the **effective stiffness matrix** $\mathbf{\hat{K}}$ and **effective load vector** $\mathbf{\hat{p}}_{n+1}$:

$$
\mathbf{\hat{K}} = \mathbf{K} + \frac{1}{\beta \Delta t^2}\mathbf{M} + \frac{\gamma}{\beta \Delta t}\mathbf{C}
$$

$$
\mathbf{\hat{p}}_{n+1} = \mathbf{p}_{n+1} + \mathbf{M} \left[ \frac{1}{\beta \Delta t^2}\mathbf{u}_n + \frac{1}{\beta \Delta t}\dot{\mathbf{u}}_n + \left(\frac{1}{2\beta} - 1\right)\ddot{\mathbf{u}}_n \right] + \mathbf{C} \left[ \frac{\gamma}{\beta \Delta t}\mathbf{u}_n + \left(\frac{\gamma}{\beta} - 1\right)\dot{\mathbf{u}}_n + \Delta t \left(\frac{\gamma}{2\beta} - 1\right)\ddot{\mathbf{u}}_n \right]
$$

The system $\mathbf{\hat{K}} \mathbf{u}_{n+1} = \mathbf{\hat{p}}_{n+1}$ is solved for $\mathbf{u}_{n+1}$, and the velocities and accelerations are subsequently updated. 

**Common Parameter Sets:**
*   **Average Acceleration:** $\beta = 1/4, \gamma = 1/2$ (Unconditionally stable, no numerical dissipation).
*   **Linear Acceleration:** $\beta = 1/6, \gamma = 1/2$ (Conditionally stable $\Delta t \le 0.551 T_n$).

## 3. Explicit Time Integration: Central Difference Method

The Central Difference method is explicit and strictly conditionally stable ($\Delta t \le \frac{T_n}{\pi}$). It requires a lumped (diagonal) mass matrix $\mathbf{M}_L$.

The finite difference approximations are:

$$
\ddot{\mathbf{u}}_n = \frac{1}{\Delta t^2} (\mathbf{u}_{n+1} - 2\mathbf{u}_n + \mathbf{u}_{n-1})
$$
$$
\dot{\mathbf{u}}_n = \frac{1}{2\Delta t} (\mathbf{u}_{n+1} - \mathbf{u}_{n-1})
$$

Yielding the explicit update equation:

$$
\left( \frac{1}{\Delta t^2}\mathbf{M}_L + \frac{1}{2\Delta t}\mathbf{C} \right) \mathbf{u}_{n+1} = \mathbf{p}_n - \left( \mathbf{K} - \frac{2}{\Delta t^2}\mathbf{M}_L \right)\mathbf{u}_n - \left( \frac{1}{\Delta t^2}\mathbf{M}_L - \frac{1}{2\Delta t}\mathbf{C} \right)\mathbf{u}_{n-1}
$$

Because the effective mass matrix is diagonal, the solution is purely vectorial (no matrix factorization required), making it extremely fast per time step for wave propagation and impact problems.

## 4. Dissipative Integration: HHT-$\alpha$ Method

The Hilber-Hughes-Taylor (HHT) method introduces numerical dissipation for high-frequency noise while retaining second-order accuracy. It modifies the equation of motion by averaging the stiffness, damping, and external force terms over the step:

$$
\mathbf{M}\ddot{\mathbf{u}}_{n+1} + (1+\alpha)\mathbf{C}\dot{\mathbf{u}}_{n+1} - \alpha\mathbf{C}\dot{\mathbf{u}}_n + (1+\alpha)\mathbf{K}\mathbf{u}_{n+1} - \alpha\mathbf{K}\mathbf{u}_n = (1+\alpha)\mathbf{p}_{n+1} - \alpha\mathbf{p}_n
$$

The Newmark parameters are strictly tied to $\alpha \in [-1/3, 0]$ to ensure unconditional stability and second-order accuracy:

$$
\beta = \frac{(1 - \alpha)^2}{4}, \qquad \gamma = \frac{1}{2} - \alpha
$$

Setting $\alpha = 0$ recovers the standard Newmark Average Acceleration method.

## 5. Frequency Response Function (FRF)

For steady-state harmonic analysis under a load $\mathbf{p}(t) = \mathbf{f} e^{i\omega t}$, the response is $\mathbf{u}(t) = \mathbf{u}_0 e^{i\omega t}$. Substituting this into the EOM yields the complex dynamic stiffness:

$$
\left(-\omega^2\mathbf{M} + i\omega\mathbf{C} + \mathbf{K}\right)\mathbf{u}_0 = \mathbf{f}
$$

The Frequency Response Function (FRF) matrix $\mathbf{H}(\omega)$ is the inverse of the dynamic stiffness:

$$
\mathbf{H}(\omega) = \left(\mathbf{K} - \omega^2\mathbf{M} + i\omega\mathbf{C}\right)^{-1}
$$

`compute_frf()` evaluates this directly over a frequency vector to construct Bode plots (magnitude and phase vs. frequency).