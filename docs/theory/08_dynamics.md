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

The consistent mass matrix is derived from the kinetic energy expression:

$$
T = \frac{1}{2} \int_{\Omega} \rho \dot{\mathbf{u}}^T \dot{\mathbf{u}} \, dV
$$

Substituting the FEM interpolation $\mathbf{u} = \mathbf{N} \mathbf{d}$:

$$
T = \frac{1}{2} \dot{\mathbf{d}}^T \underbrace{\int_{\Omega} \rho \mathbf{N}^T \mathbf{N} \, dV}_{\mathbf{M}_c} \dot{\mathbf{d}}
$$

For a T3 element, the consistent mass matrix is:

$$
\mathbf{M}_c^{T3} = \frac{\rho A}{12} \begin{bmatrix} 2 & 1 & 1 \\ 1 & 2 & 1 \\ 1 & 1 & 2 \end{bmatrix} \otimes \mathbf{I}_2
$$

For a Q4 element, numerical integration is used:

$$
\mathbf{M}_c^{Q4} = \sum_{g} w_g \rho \mathbf{N}_g^T \mathbf{N}_g |\mathbf{J}_g| t
$$

**Lumped Mass Matrix ($\mathbf{M}_L$):**

A diagonal approximation that decouples inertial forces. Two common methods:

*Row-sum lumping:*
$$
M_{L,ii} = \sum_j M_{c,ij}
$$

*HRZ (Hinton-Rock-Zienkiewicz) scaling:*
$$
M_{L, ii} = \alpha \int_{\Omega^e} \rho N_i^2 \, d\Omega \quad \text{where} \quad \alpha = \frac{\int \rho \, dV}{\sum_j \int \rho N_j^2 \, dV}
$$

Lumped mass is required for explicit time integration (Central Difference).

### 1.2. Rayleigh Damping

**Derivation:**

For mode $n$ with natural frequency $\omega_n$, the modal damping ratio is:

$$
\zeta_n = \frac{\alpha}{2\omega_n} + \frac{\beta \omega_n}{2}
$$

Given target damping ratios $\zeta_1, \zeta_2$ at frequencies $\omega_1, \omega_2$:

$$
\begin{bmatrix} \frac{1}{2\omega_1} & \frac{\omega_1}{2} \\ \frac{1}{2\omega_2} & \frac{\omega_2}{2} \end{bmatrix} \begin{bmatrix} \alpha \\ \beta \end{bmatrix} = \begin{bmatrix} \zeta_1 \\ \zeta_2 \end{bmatrix}
$$

Solving:

$$
\begin{bmatrix} \alpha \\ \beta \end{bmatrix} = \frac{2\omega_1\omega_2}{\omega_2^2 - \omega_1^2} \begin{bmatrix} \omega_2\zeta_1 - \omega_1\zeta_2 \\ \frac{\zeta_2 - \zeta_1}{\omega_1\omega_2}(\omega_2 - \omega_1) \end{bmatrix}
$$

The global damping matrix is:

$$
\mathbf{C} = \alpha \mathbf{M} + \beta \mathbf{K}
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

---

## Appendix: Mathematical Foundations

### C.1 Newmark-β derivation from Taylor series

The Newmark method can be derived by considering Taylor series expansions and introducing a weighted average of accelerations.

**Taylor expansion of displacement:**

$$
\mathbf{u}_{n+1} = \mathbf{u}_n + \Delta t \dot{\mathbf{u}}_n + \frac{\Delta t^2}{2} \ddot{\mathbf{u}}_n + \frac{\Delta t^3}{6} \dddot{\mathbf{u}}_n + O(\Delta t^4)
$$

**Newmark assumption:** Replace the unknown third derivative term with a weighted average:

$$
\frac{\Delta t^2}{2} \ddot{\mathbf{u}}_n + \frac{\Delta t^3}{6} \dddot{\mathbf{u}}_n \approx \frac{\Delta t^2}{2} \left[ (1-2\beta)\ddot{\mathbf{u}}_n + 2\beta \ddot{\mathbf{u}}_{n+1} \right]
$$

This gives the Newmark displacement equation. The velocity equation follows similarly:

$$
\dot{\mathbf{u}}_{n+1} = \dot{\mathbf{u}}_n + \Delta t \ddot{\mathbf{u}}_n + \frac{\Delta t^2}{2} \dddot{\mathbf{u}}_n + O(\Delta t^3)
$$

Approximating:

$$
\ddot{\mathbf{u}}_n + \frac{\Delta t}{2} \dddot{\mathbf{u}}_n \approx (1-\gamma)\ddot{\mathbf{u}}_n + \gamma \ddot{\mathbf{u}}_{n+1}
$$

### C.2 Stability analysis via amplification matrix

The stability of Newmark integration is analyzed through the amplification matrix. For a single DOF system:

$$
m\ddot{u} + c\dot{u} + ku = p(t)
$$

Define the state vector $\mathbf{x} = [u, \dot{u}]^T$. The recurrence relation is:

$$
\mathbf{x}_{n+1} = \mathbf{A} \mathbf{x}_n + \mathbf{L} p_{n+1}
$$

**Spectral radius condition:** For unconditional stability:

$$
\rho(\mathbf{A}) \leq 1 \quad \forall \, \Omega = \omega \Delta t
$$

**Proof for average acceleration ($\beta = 1/4$, $\gamma = 1/2$):**

The amplification matrix for the undamped case is:

$$
\mathbf{A} = \begin{bmatrix} 
\frac{1 - \Omega^2/4}{1 + \Omega^2/4} & \frac{\Delta t}{1 + \Omega^2/4} \\
-\frac{\Omega^2}{\Delta t(1 + \Omega^2/4)} & \frac{1 - \Omega^2/4}{1 + \Omega^2/4}
\end{bmatrix}
$$

The eigenvalues are:

$$
\lambda_{1,2} = \frac{(1 - \Omega^2/4) \pm i\Omega}{1 + \Omega^2/4}
$$

The spectral radius:

$$
\rho = |\lambda| = \sqrt{\frac{(1 - \Omega^2/4)^2 + \Omega^2}{(1 + \Omega^2/4)^2}} = 1
$$

for all $\Omega > 0$. This proves unconditional stability with no amplitude decay. ∎

### C.3 Accuracy order analysis

**Local truncation error:** For the Newmark displacement equation:

$$
\tau = \mathbf{u}(t_{n+1}) - \mathbf{u}_n - \Delta t \dot{\mathbf{u}}_n - \frac{\Delta t^2}{2}\left[(1-2\beta)\ddot{\mathbf{u}}_n + 2\beta \ddot{\mathbf{u}}_{n+1}\right]
$$

Expanding $\ddot{\mathbf{u}}_{n+1}$ in Taylor series:

$$
\ddot{\mathbf{u}}_{n+1} = \ddot{\mathbf{u}}_n + \Delta t \dddot{\mathbf{u}}_n + O(\Delta t^2)
$$

Substituting:

$$
\tau = \frac{\Delta t^3}{6}\dddot{\mathbf{u}}_n - \beta \Delta t^3 \dddot{\mathbf{u}}_n + O(\Delta t^4) = \Delta t^3 \left(\frac{1}{6} - \beta\right)\dddot{\mathbf{u}}_n + O(\Delta t^4)
$$

For $\gamma = 1/2$, the method is **second-order accurate** regardless of $\beta$.

For $\gamma \neq 1/2$, artificial damping is introduced:

$$
\zeta_{num} = \frac{\gamma - 1/2}{2} \omega \Delta t
$$

### C.4 Critical time step for explicit methods

For the central difference method, the critical time step is:

$$
\Delta t_{cr} = \frac{2}{\omega_{max}}
$$

where $\omega_{max}$ is the highest natural frequency of the system.

**Derivation:** For the undamped recurrence:

$$
u_{n+1} = \left(2 - \Omega^2\right)u_n - u_{n-1}
$$

The characteristic equation is:

$$
\lambda^2 - (2 - \Omega^2)\lambda + 1 = 0
$$

For stability, we need $|\lambda| \leq 1$, which requires:

$$
|2 - \Omega^2| \leq 2 \implies 0 \leq \Omega^2 \leq 4 \implies \omega \Delta t \leq 2
$$

Therefore: $\Delta t \leq \frac{2}{\omega_{max}}$ ∎

**Practical estimate:** For FEM meshes:

$$
\omega_{max} \approx \frac{2c}{h_{min}}
$$

where $c = \sqrt{E/\rho}$ is the wave speed and $h_{min}$ is the smallest element dimension.

Thus:

$$
\Delta t_{cr} \approx \frac{h_{min}}{c}
$$

### C.5 Modal superposition

For large systems, direct integration is expensive. Modal superposition reduces the problem size.

**Transformation:** Using mass-normalized modes $\boldsymbol{\Phi}$:

$$
\mathbf{u} = \boldsymbol{\Phi} \mathbf{q}
$$

The equations decouple:

$$
\ddot{q}_n + 2\zeta_n \omega_n \dot{q}_n + \omega_n^2 q_n = \phi_n^T \mathbf{p}(t)
$$

Each modal equation is a SDOF system that can be solved analytically or numerically.

**Truncation:** Only modes with $\omega_n \leq \omega_{cut}$ are retained. The error is bounded by:

$$
\|\mathbf{u} - \mathbf{u}_{approx}\| \leq \frac{\|\mathbf{p}\|}{\omega_{cut}^2 m_{eff}}
$$

where $m_{eff}$ is the effective modal mass of neglected modes.

### C.6 Energy conservation

For conservative systems ($\mathbf{C} = 0$, $\mathbf{p} = 0$), the total energy should be preserved.

**Definition:**

$$
E_n = \frac{1}{2}\dot{\mathbf{u}}_n^T \mathbf{M} \dot{\mathbf{u}}_n + \frac{1}{2}\mathbf{u}_n^T \mathbf{K} \mathbf{u}_n
$$

**Energy behavior:**

| Method | Energy |
|--------|--------|
| Average acceleration | Conserved exactly |
| Linear acceleration | Grows (unstable for large $\Delta t$) |
| Central difference | Bounded oscillation |
| HHT-α | Dissipated (by design) |

The average acceleration method satisfies:

$$
E_{n+1} = E_n \quad \forall n
$$

This is why it's preferred for long-duration simulations where energy drift is unacceptable.