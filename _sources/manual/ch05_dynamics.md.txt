# Chapter 5: Structural Dynamics

Structural dynamics extends standard statics by incorporating inertial (mass) and damping forces into the equations of motion:
$$ \mathbf{M} \ddot{\mathbf{u}}(t) + \mathbf{C} \dot{\mathbf{u}}(t) + \mathbf{K} \mathbf{u}(t) = \mathbf{P}(t) $$

## 5.1 Mass Matrices

The global mass matrix $\mathbf{M}$ represents the inertial distribution. `femlabpy` provides two formulations:

1. **Consistent Mass Matrix:** Derived using the same shape functions $\mathbf{N}$ as the stiffness matrix.
   $$ \mathbf{M}_e^c = \int_V \rho \mathbf{N}^T \mathbf{N} \, dV $$
   This creates a coupled, banded matrix that perfectly preserves the total kinetic energy representation.

2. **Lumped Mass Matrix:** A strictly diagonal matrix obtained by summing the rows of the consistent mass matrix (Row-Sum technique) or using specialized numerical integration. Diagonal mass matrices are computationally cheap and strictly required for explicit time-integration algorithms.

## 5.2 Modal Analysis

Free, undamped vibration governed by $\mathbf{K} \mathbf{\phi} = \omega^2 \mathbf{M} \mathbf{\phi}$ is solved via `solve_modal()`.

- **Natural Frequencies ($f_n$):** Extracted from the eigenvalues $\omega_n = 2\pi f_n$.
- **Mode Shapes ($\phi_n$):** The eigenvectors, mass-normalized such that $\phi_n^T \mathbf{M} \phi_n = 1$.

### Effective Modal Mass
To understand how much of the total mass is activated by a specific mode $n$ in a spatial direction $j$, we compute the participation factor $\Gamma_{nj}$ and the effective modal mass $m_{eff}$:
$$ \Gamma_{nj} = \frac{\phi_n^T \mathbf{M} \mathbf{r}_j}{\phi_n^T \mathbf{M} \phi_n} $$
$$ m_{eff,nj} = \Gamma_{nj}^2 \times (\phi_n^T \mathbf{M} \phi_n) $$
where $\mathbf{r}_j$ is the influence vector for direction $j$.

## 5.3 Damping Models

Damping $\mathbf{C}$ is rarely derived explicitly from element geometry; instead, it is formulated at the global level.

### Rayleigh Damping
The most common approach is proportional damping:
$$ \mathbf{C} = \alpha \mathbf{M} + \beta \mathbf{K} $$
Given two critical damping ratios $\zeta_1, \zeta_2$ at circular frequencies $\omega_1, \omega_2$, the multipliers are found by solving:
$$ \zeta_n = \frac{\alpha}{2\omega_n} + \frac{\beta \omega_n}{2} $$
The `rayleigh_coefficients` and `rayleigh_damping` functions automate this.

## 5.4 Time Integration

To solve the differential equations under arbitrary loads $\mathbf{P}(t)$, we use step-by-step numerical integration.

### Newmark-$\beta$ Method (Implicit)
The unconditionally stable average acceleration method assumes acceleration is constant over a time step $\Delta t$ ($\gamma=0.5, \beta=0.25$).
The displacement and velocity updates are:
$$ \mathbf{u}_{t+\Delta t} = \mathbf{u}_t + \Delta t \dot{\mathbf{u}}_t + \Delta t^2 \left( (0.5 - \beta)\ddot{\mathbf{u}}_t + \beta \ddot{\mathbf{u}}_{t+\Delta t} \right) $$
$$ \dot{\mathbf{u}}_{t+\Delta t} = \dot{\mathbf{u}}_t + \Delta t \left( (1-\gamma)\ddot{\mathbf{u}}_t + \gamma \ddot{\mathbf{u}}_{t+\Delta t} \right) $$
At each step, an effective stiffness matrix $\mathbf{K}_{eff}$ is factored and solved:
$$ \mathbf{K}_{eff} = \mathbf{K} + a_0 \mathbf{M} + a_1 \mathbf{C} $$
where $a_0 = \frac{1}{\beta \Delta t^2}$ and $a_1 = \frac{\gamma}{\beta \Delta t}$.

### Central Difference Method (Explicit)
For impact or wave propagation, the explicit central difference method evaluates equilibrium at time $t$ without factoring a global matrix, provided $\mathbf{M}$ is diagonal. However, it is conditionally stable; the time step must strictly obey $\Delta t \le \Delta t_{critical} = \frac{2}{\omega_{max}}$.
