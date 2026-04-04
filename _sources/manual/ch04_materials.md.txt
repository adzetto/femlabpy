# Chapter 4: Material Models

`femlabpy` supports a variety of constitutive models, mapped through the material matrix $\mathbf{D}$ relating stress to strain: $\sigma = \mathbf{D} \epsilon$. This chapter outlines the theoretical foundations for elasticity and plasticity implementations.

## 4.1 Linear Elasticity

For an isotropic linear elastic material, the compliance tensor depends solely on Young's Modulus $E$ and Poisson's ratio $\nu$. In 2D formulations, we distinguish between two states:

### Plane Stress
Assume $\sigma_z = \tau_{xz} = \tau_{yz} = 0$. The material stiffness matrix is:
$$ \mathbf{D} = \frac{E}{1 - \nu^2} \begin{bmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \frac{1-\nu}{2}
\end{bmatrix} $$

### Plane Strain
Assume $\epsilon_z = \gamma_{xz} = \gamma_{yz} = 0$. The material stiffness matrix is:
$$ \mathbf{D} = \frac{E}{(1+\nu)(1-2\nu)} \begin{bmatrix}
1-\nu & \nu & 0 \\
\nu & 1-\nu & 0 \\
0 & 0 & \frac{1-2\nu}{2}
\end{bmatrix} $$

## 4.2 von Mises Elastoplasticity

For ductile metals, the von Mises yield criterion is employed with isotropic hardening.

### Yield Function
The yield function $f$ defines the elastic limit based on the deviatoric stress invariants:
$$ f(\sigma, \bar{\epsilon}_p) = \sigma_{eq} - (\sigma_y + H \bar{\epsilon}_p) \le 0 $$
where $\sigma_{eq} = \sqrt{\frac{3}{2} \mathbf{s}:\mathbf{s}}$ is the equivalent von Mises stress, $\mathbf{s}$ is the deviatoric stress tensor, $\sigma_y$ is the initial yield stress, $H$ is the hardening modulus, and $\bar{\epsilon}_p$ is the accumulated equivalent plastic strain.

### Radial Return Mapping Algorithm
When the elastic trial stress violates the yield condition ($f_{trial} > 0$), a backward-Euler integration scheme (Radial Return) is executed at the Gauss point level to return the stress to the yield surface:

1. Compute elastic trial stress: $\sigma_{trial} = \mathbf{D} \epsilon_{total}$
2. Compute trial deviatoric stress $\mathbf{s}_{trial}$ and equivalent stress $q_{trial}$.
3. Check yield condition. If plastic, compute the plastic multiplier increment $\Delta \gamma$:
   $$ \Delta \gamma = \frac{q_{trial} - \sigma_y - H \bar{\epsilon}_p}{3G + H} $$
   where $G = \frac{E}{2(1+\nu)}$ is the shear modulus.
4. Update the stress tensor:
   $$ \sigma_{n+1} = \mathbf{s}_{trial} \left( 1 - \frac{3G \Delta \gamma}{q_{trial}} \right) + p \mathbf{I} $$
   where $p$ is the hydrostatic pressure.
5. Update plastic strain: $\bar{\epsilon}_{p, n+1} = \bar{\epsilon}_{p, n} + \Delta \gamma$.

### Consistent Tangent Modulus
To preserve the quadratic convergence of the global Newton-Raphson iterations, the algorithmic consistent tangent stiffness $\mathbf{D}^{ep}$ must be computed, modifying the elastic $\mathbf{D}$ matrix based on the normal flow direction $\mathbf{n} = \frac{3}{2} \frac{\mathbf{s}}{q}$.

## 4.3 Drucker-Prager Plasticity

For soils and concrete, the yield strength depends heavily on hydrostatic pressure. The Drucker-Prager model approximates the Mohr-Coulomb failure surface using the invariants $I_1$ (pressure) and $J_2$ (deviatoric shear).

$$ f(\sigma) = \sqrt{J_2} + \alpha I_1 - k \le 0 $$

The parameters $\alpha$ and $k$ are derived from the internal friction angle $\phi$ and cohesion $c$. The return mapping algorithm iteratively projects the trial stress back to the conical yield surface or to the cone apex if tension prevails.
