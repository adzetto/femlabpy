# Chapter 8: Custom Element Development

The architecture of `femlabpy` makes it exceptionally straightforward to write and integrate your own finite elements. Unlike object-oriented frameworks where you must subclass deeply nested structures, `femlabpy` requires exactly two mathematically pure functions per element type:

1. **Stiffness Routine (`ke...`)**: A function that returns the $n \times n$ elemental stiffness matrix.
2. **Force Routine (`qe...`)**: A function that returns the $n \times 1$ elemental internal force vector (and optionally, stresses and strains).

This chapter walks you through the complete implementation of a **2D Euler-Bernoulli Beam Element** from scratch, demonstrating how it hooks natively into the existing global assembly engines (`assmk` and `assmq`).

## 8.1 Theoretical Formulation (2D Beam)

A 2D Euler-Bernoulli beam element has two nodes, with 3 Degrees of Freedom (DOFs) per node:
- $u_x$: Axial displacement
- $u_y$: Transverse deflection
- $\theta_z$: Rotation

### Stiffness Matrix
The local $6 \times 6$ stiffness matrix for a beam of length $L$, Young's Modulus $E$, Area $A$, and Moment of Inertia $I$ is well-known:

$$
\mathbf{K}_{local} = \begin{bmatrix}
\frac{EA}{L} & 0 & 0 & -\frac{EA}{L} & 0 & 0 \\
0 & \frac{12EI}{L^3} & \frac{6EI}{L^2} & 0 & -\frac{12EI}{L^3} & \frac{6EI}{L^2} \\
0 & \frac{6EI}{L^2} & \frac{4EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{2EI}{L} \\
-\frac{EA}{L} & 0 & 0 & \frac{EA}{L} & 0 & 0 \\
0 & -\frac{12EI}{L^3} & -\frac{6EI}{L^2} & 0 & \frac{12EI}{L^3} & -\frac{6EI}{L^2} \\
0 & \frac{6EI}{L^2} & \frac{2EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{4EI}{L}
\end{bmatrix}
$$

### Coordinate Transformation
To orient the beam in 2D space at an angle $\alpha$, we apply the $6 \times 6$ transformation matrix $\mathbf{R}$:
$$ \mathbf{K}_{global} = \mathbf{R}^T \mathbf{K}_{local} \mathbf{R} $$
$$ \mathbf{R} = \begin{bmatrix} \cos\alpha & \sin\alpha & 0 & 0 & 0 & 0 \\ -\sin\alpha & \cos\alpha & 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 & 0 \\ \vdots & & & \ddots & & \vdots \end{bmatrix} $$

## 8.2 Implementing the Core Routines

### 1. Element Stiffness Function (`kebeam2d`)

```python
import numpy as np

def kebeam2d(Xe, Ge):
    """
    Compute the 6x6 stiffness matrix for a 2D beam element.

    Parameters
    ----------
    Xe : ndarray
        2x2 array of node coordinates: [[x1, y1], [x2, y2]]
    Ge : array_like
        Material properties: [Area, Moment of Inertia (I), Young's Modulus (E)]

    Returns
    -------
    Ke : ndarray
        6x6 global element stiffness matrix.
    """
    A, I, E = Ge[0], Ge[1], Ge[2]
    
    # Calculate length and orientation
    dx = Xe[1, 0] - Xe[0, 0]
    dy = Xe[1, 1] - Xe[0, 1]
    L = np.sqrt(dx**2 + dy**2)
    
    c = dx / L
    s = dy / L
    
    # Local stiffness matrix
    k_axial = E * A / L
    k_shear = 12 * E * I / L**3
    k_mom1  = 6 * E * I / L**2
    k_mom2  = 4 * E * I / L
    k_mom3  = 2 * E * I / L
    
    K_loc = np.array([
        [ k_axial,  0,        0,       -k_axial,  0,        0      ],
        [ 0,        k_shear,  k_mom1,   0,       -k_shear,  k_mom1 ],
        [ 0,        k_mom1,   k_mom2,   0,       -k_mom1,   k_mom3 ],
        [-k_axial,  0,        0,        k_axial,  0,        0      ],
        [ 0,       -k_shear, -k_mom1,   0,        k_shear, -k_mom1 ],
        [ 0,        k_mom1,   k_mom3,   0,       -k_mom1,   k_mom2 ]
    ])
    
    # Transformation matrix
    R = np.array([
        [ c,  s,  0,  0,  0,  0],
        [-s,  c,  0,  0,  0,  0],
        [ 0,  0,  1,  0,  0,  0],
        [ 0,  0,  0,  c,  s,  0],
        [ 0,  0,  0, -s,  c,  0],
        [ 0,  0,  0,  0,  0,  1]
    ])
    
    # K_global = R^T * K_local * R
    Ke = R.T @ K_loc @ R
    return Ke
```

### 2. Element Internal Force Function (`qebeam2d`)

```python
def qebeam2d(Xe, Ge, Ue):
    """
    Compute internal forces for a 2D beam element.

    Parameters
    ----------
    Ue : ndarray
        6x1 vector of nodal displacements in global coordinates.

    Returns
    -------
    qe : ndarray
        6x1 internal force vector in global coordinates.
    N, V, M : floats
        Internal axial force, shear force, and bending moment.
    """
    Ke = kebeam2d(Xe, Ge)
    
    # Internal force in global coordinates
    qe = Ke @ Ue.reshape(6, 1)
    
    # Transform displacements to local to get intuitive engineering forces
    # ... (transformation omitted for brevity)
    
    return qe, 0.0, 0.0, 0.0  # Simplification for tutorial
```

## 8.3 The Global Driver: Hooking into `assmk`

To assemble a mesh of 1,000 beams, you don't write `for` loops in your main script. You create a single "driver" function, `kbeam2d`, which utilizes the extremely fast `femlabpy.assmk` assembly engine.

```python
from femlabpy.assembly import assmk
from femlabpy._helpers import is_sparse

def kbeam2d(K, T, X, G):
    """
    Assemble the global stiffness matrix for all 2D beam elements.

    Parameters
    ----------
    K : ndarray or scipy.sparse matrix
        The un-assembled or partially assembled global stiffness matrix.
    T : ndarray
        Topology matrix: [node1, node2, prop_id]
    X : ndarray
        Node coordinates: [x, y]
    G : ndarray
        Material property table.

    Returns
    -------
    K : Assembled global stiffness matrix.
    """
    nel = T.shape[0]
    dof = 3  # 3 DOFs per node (ux, uy, theta_z)
    
    for e in range(nel):
        # 1-based node IDs from the topology table
        n1, n2 = int(T[e, 0]), int(T[e, 1])
        prop_id = int(T[e, 2]) - 1
        
        # Extract coordinates and properties
        Xe = np.vstack((X[n1-1], X[n2-1]))
        Ge = G[prop_id]
        
        # Compute 6x6 element matrix
        Ke = kebeam2d(Xe, Ge)
        
        # Assemble using the native femlabpy engine
        nodes = [n1, n2]
        if is_sparse(K):
            # Specialized fast sparse assembly (internal magic)
            # K is modified in-place if using LIL format
            K = assmk(K, Ke, nodes, dof)
        else:
            K = assmk(K, Ke, nodes, dof)
            
    return K
```

That's it! Your custom element is now a first-class citizen in the `femlabpy` ecosystem. You can mix it with standard 2D Quads, apply `setbc` constraints to the 3rd DOF (rotations), and solve dynamic matrices using `solve_newmark`. The open, mathematically-driven architecture makes extending the library trivial.