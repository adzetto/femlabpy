# Chapter 1: Fundamentals & Data Structures

## 1.1 The femlabpy Philosophy

`femlabpy` is a pure-Python finite element library that inherits the array-based, vectorized philosophy of the MATLAB/Scilab FemLab toolboxes. Unlike many modern object-oriented finite element frameworks that represent nodes, elements, and materials as deeply nested class instances, `femlabpy` relies entirely on flat, contiguous `numpy` arrays. This design choice has several profound implications:

1. **Performance:** By storing data in NumPy arrays, we leverage highly optimized C and Fortran backends for matrix operations, assembly, and slicing.
2. **Transparency:** The mathematical mapping from theory to code is direct. A stiffness matrix $K_e$ is simply a 2D array; a displacement vector $u$ is a 1D array.
3. **Didactic Clarity:** Students and researchers can inspect the entire state of the model at any point by simply printing an array. There are no hidden states or complex getters/setters.

## 1.2 Core Data Structures

A finite element model in `femlabpy` is completely described by five core matrices. Understanding the mathematical and structural shape of these arrays is the fundamental prerequisite for using the library.

### 1.2.1 Node Coordinates Matrix (`X`)

The matrix `X` stores the spatial coordinates of all nodes in the global coordinate system. For a 2D problem with $N$ nodes, `X` is an $N \times 2$ matrix. The row index corresponds to the global node ID.

$$
\mathbf{X} = \begin{bmatrix}
x_1 & y_1 \\
x_2 & y_2 \\
\vdots & \vdots \\
x_N & y_N
\end{bmatrix}
$$

*Note: `femlabpy` strictly uses 1-based indexing for node references in topology arrays to maintain compatibility with legacy meshes and mathematical conventions, although the Python arrays themselves are 0-indexed under the hood.*

### 1.2.2 Topology Matrix (`T`)

The topology matrix `T` defines the connectivity of the elements. For a mesh of 4-node quadrilaterals (Q4) with $E$ elements, `T` is an $E \times 5$ integer array. The first 4 columns contain the 1-based global node IDs that make up the element, ordered counter-clockwise. The final column is the material property ID.

$$
\mathbf{T} = \begin{bmatrix}
n_{1,1} & n_{1,2} & n_{1,3} & n_{1,4} & \text{prop}_1 \\
\vdots & \vdots & \vdots & \vdots & \vdots \\
n_{E,1} & n_{E,2} & n_{E,3} & n_{E,4} & \text{prop}_E
\end{bmatrix}
$$

### 1.2.3 Material Properties (`G`)

The material matrix `G` contains the physical properties assigned to the elements. For 2D plane stress/strain problems, a row typically takes the form: `[E, \nu, \text{type}, t, \rho]`, where:
- $E$: Young's Modulus
- $\nu$: Poisson's Ratio
- $\text{type}$: 1 for Plane Stress, 2 for Plane Strain
- $t$: Thickness (for Plane Stress)
- $\rho$: Density (for dynamic mass matrices)

### 1.2.4 Boundary Constraints (`C`)

Dirichlet boundary conditions (prescribed displacements) are stored in the `C` array. Each row specifies a constrained node, the degree of freedom (DOF) index (1 for $u_x$, 2 for $u_y$), and the prescribed value $\bar{u}$.

$$
\mathbf{C} = \begin{bmatrix}
\text{node\_id}_1 & \text{dof}_1 & \bar{u}_1 \\
\vdots & \vdots & \vdots
\end{bmatrix}
$$

### 1.2.5 Point Loads (`P`)

Neumann boundary conditions (point loads) are stored in the `P` array, following a similar format to `C`: `[node_id, dof, force_value]`.

## 1.3 The Finite Element Analysis Sequence

The standard linear static analysis in `femlabpy` follows a rigorous mathematical sequence:

1. **Initialization:** The global stiffness matrix $\mathbf{K}$ and load vector $\mathbf{p}$ are initialized to zero.
   ```python
   K, p = fp.init(nn, dof)
   ```
2. **Assembly:** Element stiffness matrices $\mathbf{K}_e$ are computed and assembled into $\mathbf{K}$.
   $$ \mathbf{K} = \sum_{e=1}^{E} \mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e $$
   ```python
   K = fp.kq4e(K, T, X, G)
   ```
3. **Load Application:** Point loads are mapped into the global load vector $\mathbf{p}$.
   ```python
   p = fp.setload(p, P, dof)
   ```
4. **Boundary Conditions:** Constraints are applied using the penalty method, modifying $\mathbf{K}$ and $\mathbf{p}$ in place.
   ```python
   K_bc, p_bc, _ = fp.setbc(K, p, C, dof)
   ```
5. **Solution:** The algebraic system $\mathbf{K}_{bc} \mathbf{u} = \mathbf{p}_{bc}$ is solved for the nodal displacements $\mathbf{u}$.
   ```python
   u = np.linalg.solve(K_bc, p_bc)
   ```
6. **Internal Forces Recovery:** Stresses and strains are computed at the element Gauss points.
   ```python
   q, S, E = fp.qq4e(np.zeros_like(p), T, X, G, u)
   ```