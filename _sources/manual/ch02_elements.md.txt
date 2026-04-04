# Chapter 2: Element Library

The element library in `femlabpy` provides the mathematical mapping from the continuous differential equations of solid mechanics to the discrete algebraic space. This chapter details the shape functions, integration schemes, and strain-displacement matrices ($\mathbf{B}$) used for our core elements.

## 2.1 1D Bar Element (`bar`)

The bar element supports axial tension and compression. In `femlabpy`, it is formulated to handle large deformations by utilizing the Green-Lagrange strain measure.

### Kinematics and Strain
For large deformations, the axial strain $\epsilon$ is defined as:
$$ \epsilon = \frac{du}{dx} + \frac{1}{2}\left(\frac{du}{dx}\right)^2 $$

The internal virtual work gives rise to two stiffness components: the standard linear elastic material stiffness $\mathbf{K}_m$ and the geometric stiffness $\mathbf{K}_g$, which accounts for the effect of the internal axial force $N$ on the transverse stiffness.

$$ \mathbf{K}_e = \mathbf{K}_m + \mathbf{K}_g $$
$$ \mathbf{K}_m = \frac{EA}{L} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}, \quad \mathbf{K}_g = \frac{N}{L} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix} $$

## 2.2 Constant Strain Triangle (`t3`)

The 3-node triangular element (CST) is the simplest 2D element. The displacement field is interpolated linearly using area coordinates $L_1, L_2, L_3$.

$$ u(x,y) = N_1 u_1 + N_2 u_2 + N_3 u_3 $$

Since the shape functions $N_i$ are linear, their spatial derivatives are constant over the element. Consequently, the strain-displacement matrix $\mathbf{B}$ is constant:

$$ \mathbf{B} = \frac{1}{2A} \begin{bmatrix}
y_{23} & 0 & y_{31} & 0 & y_{12} & 0 \\
0 & x_{32} & 0 & x_{13} & 0 & x_{21} \\
x_{32} & y_{23} & x_{13} & y_{31} & x_{21} & y_{12}
\end{bmatrix} $$

where $x_{ij} = x_i - x_j$ and $y_{ij} = y_i - y_j$, and $A$ is the element area. The element stiffness matrix is derived exactly without numerical integration:
$$ \mathbf{K}_e = \mathbf{B}^T \mathbf{D} \mathbf{B} A t $$

## 2.3 Isoparametric Quadrilateral (`q4`)

The 4-node quadrilateral element uses a bilinear isoparametric formulation. The element geometry and displacements are mapped from a natural coordinate system $(\xi, \eta) \in [-1, 1]$ to the physical system $(x, y)$.

### Shape Functions
The bilinear shape functions are:
$$ N_i(\xi, \eta) = \frac{1}{4} (1 + \xi_i \xi) (1 + \eta_i \eta) $$

### The Jacobian Matrix
To compute spatial derivatives with respect to $x$ and $y$, we apply the chain rule via the Jacobian matrix $\mathbf{J}$:

$$ \mathbf{J} = \begin{bmatrix}
\frac{\partial x}{\partial \xi} & \frac{\partial y}{\partial \xi} \\
\frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta}
\end{bmatrix} = \sum_{i=1}^4 \begin{bmatrix}
\frac{\partial N_i}{\partial \xi} x_i & \frac{\partial N_i}{\partial \xi} y_i \\
\frac{\partial N_i}{\partial \eta} x_i & \frac{\partial N_i}{\partial \eta} y_i
\end{bmatrix} $$

The Cartesian shape function derivatives are then:
$$ \begin{Bmatrix} \frac{\partial N_i}{\partial x} \\ \frac{\partial N_i}{\partial y} \end{Bmatrix} = \mathbf{J}^{-1} \begin{Bmatrix} \frac{\partial N_i}{\partial \xi} \\ \frac{\partial N_i}{\partial \eta} \end{Bmatrix} $$

### Gauss Quadrature
The stiffness matrix requires integrating over the element area. We use $2 \times 2$ Gauss-Legendre quadrature:

$$ \mathbf{K}_e = \int_{-1}^{1} \int_{-1}^{1} \mathbf{B}^T \mathbf{D} \mathbf{B} |\mathbf{J}| t \, d\xi \, d\eta \approx \sum_{i=1}^{2} \sum_{j=1}^{2} w_i w_j \mathbf{B}^T(\xi_i, \eta_j) \mathbf{D} \mathbf{B}(\xi_i, \eta_j) |\mathbf{J}(\xi_i, \eta_j)| t $$

where the integration points are $\xi_i, \eta_j \in \pm \frac{1}{\sqrt{3}}$ and weights $w_i = 1.0$.

## 2.4 3D Elements (`t4`, `h8`)

The concepts extend naturally to three dimensions. The 8-node hexahedron (`h8`) utilizes trilinear shape functions $N_i(\xi, \eta, \zeta)$ and employs $2 \times 2 \times 2$ Gauss integration. The Jacobian becomes a $3 \times 3$ matrix, and the strain-displacement matrix $\mathbf{B}$ expands to compute the 6 components of the 3D strain tensor.
