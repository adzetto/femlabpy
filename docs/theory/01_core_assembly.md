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

# Core Assembly in Finite Element Analysis

## Introduction

Finite element assembly is the step that turns element-level quantities into the global linear system solved by `femlabpy`. This chapter provides the rigorous mathematical foundation for the assembly process, deriving it from variational principles.

---

## 1. Variational Foundation

### 1.1 Principle of minimum potential energy

The finite element method is rooted in the principle of minimum potential energy. For an elastic body, the total potential energy functional is:

$$
\Pi(\mathbf{u}) = U(\mathbf{u}) - W(\mathbf{u})
$$

where:
- $U(\mathbf{u}) = \frac{1}{2}\int_\Omega \boldsymbol{\varepsilon}^T \mathbf{D} \boldsymbol{\varepsilon} \, dV$ is the strain energy
- $W(\mathbf{u}) = \int_\Omega \mathbf{u}^T \mathbf{b} \, dV + \int_{\Gamma_t} \mathbf{u}^T \mathbf{t} \, dS$ is the work of external forces

The equilibrium configuration minimizes $\Pi$:

$$
\delta \Pi = 0 \quad \forall \, \delta\mathbf{u} \text{ admissible}
$$

### 1.2 Weak form derivation

Taking the first variation of the potential energy:

$$
\delta\Pi = \int_\Omega \delta\boldsymbol{\varepsilon}^T \boldsymbol{\sigma} \, dV - \int_\Omega \delta\mathbf{u}^T \mathbf{b} \, dV - \int_{\Gamma_t} \delta\mathbf{u}^T \mathbf{t} \, dS = 0
$$

This is the **weak form** or **principle of virtual work**:

$$
\int_\Omega \delta\boldsymbol{\varepsilon}^T \boldsymbol{\sigma} \, dV = \int_\Omega \delta\mathbf{u}^T \mathbf{b} \, dV + \int_{\Gamma_t} \delta\mathbf{u}^T \mathbf{t} \, dS
$$

**Proof of equivalence to strong form:**

Using the relation $\delta\boldsymbol{\varepsilon} = \mathbf{L}\delta\mathbf{u}$ where $\mathbf{L}$ is the differential operator, and integrating by parts:

$$
\int_\Omega \delta\mathbf{u}^T \mathbf{L}^T \boldsymbol{\sigma} \, dV = -\int_\Omega \delta\mathbf{u}^T (\nabla \cdot \boldsymbol{\sigma}) \, dV + \int_{\Gamma} \delta\mathbf{u}^T (\boldsymbol{\sigma} \cdot \mathbf{n}) \, dS
$$

Since $\delta\mathbf{u}$ is arbitrary, this implies:
- $\nabla \cdot \boldsymbol{\sigma} + \mathbf{b} = \mathbf{0}$ in $\Omega$ (equilibrium)
- $\boldsymbol{\sigma} \cdot \mathbf{n} = \mathbf{t}$ on $\Gamma_t$ (traction BC)

### 1.3 Galerkin discretization

We approximate the displacement field using shape functions:

$$
\mathbf{u}(\mathbf{x}) \approx \mathbf{u}^h(\mathbf{x}) = \sum_{i=1}^{n} N_i(\mathbf{x}) \mathbf{d}_i = \mathbf{N}(\mathbf{x}) \mathbf{d}
$$

The Galerkin method requires:

$$
\delta\mathbf{u}^h \in \text{span}\{N_1, N_2, \ldots, N_n\}
$$

Substituting into the weak form and using $\delta\boldsymbol{\varepsilon} = \mathbf{B}\delta\mathbf{d}$:

$$
\delta\mathbf{d}^T \left[ \int_\Omega \mathbf{B}^T \mathbf{D} \mathbf{B} \, dV \right] \mathbf{d} = \delta\mathbf{d}^T \left[ \int_\Omega \mathbf{N}^T \mathbf{b} \, dV + \int_{\Gamma_t} \mathbf{N}^T \mathbf{t} \, dS \right]
$$

Since this must hold for all $\delta\mathbf{d}$:

$$
\mathbf{K}\mathbf{d} = \mathbf{F}
$$

where $\mathbf{K} = \int_\Omega \mathbf{B}^T \mathbf{D} \mathbf{B} \, dV$ and $\mathbf{F}$ is the load vector.

---

## 2. Assembly Identity and Derivation

The global stiffness matrix $\mathbf{K}$ and global internal force vector $\mathbf{q}$ are constructed by summing the contributions from all $N_e$ elements in the mesh. Let the global displacement vector be $\mathbf{U}$ (size $N \times 1$) and the local displacement vector for element $e$ be $\mathbf{u}_e$ (size $n_e \times 1$).

### Boolean connectivity matrix

The relationship between the global degrees of freedom and the local degrees of freedom for an element $e$ can be written with a Boolean connectivity matrix $\mathbf{L}_e$. Its role is to pick the global entries that belong to the element.

The gather operation for local displacements $\mathbf{u}_e$ is:

$$
\mathbf{u}_e = \mathbf{L}_e \mathbf{U}
$$

where the entries of $\mathbf{L}_e$ are defined as:

$$
(\mathbf{L}_e)_{ij} = 
\begin{cases} 
1 & \text{if local DOF } i \text{ corresponds to global DOF } j \\
0 & \text{otherwise}
\end{cases}
$$

Because each local DOF maps to exactly one global DOF, each row of $\mathbf{L}_e$ contains exactly one `1` and all other entries are `0`.

The transpose $\mathbf{L}_e^T$ performs the reverse operation. It scatters element contributions back into the global array. That is the algebraic origin of the assembly operation used in software.

### Strain energy and the assembly operator

To derive the global assembly equation for the stiffness matrix, consider the total strain energy of the system, $U_{\text{total}}$, which is the sum of the strain energies of the individual elements:

$$
U_{\text{total}} = \sum_{e=1}^{N_e} U_e = \sum_{e=1}^{N_e} \frac{1}{2} \mathbf{u}_e^T \mathbf{K}_e \mathbf{u}_e
$$

Substitute the Boolean mapping $\mathbf{u}_e = \mathbf{L}_e \mathbf{U}$ into the energy equation:

$$
U_{\text{total}} = \sum_{e=1}^{N_e} \frac{1}{2} (\mathbf{L}_e \mathbf{U})^T \mathbf{K}_e (\mathbf{L}_e \mathbf{U}) = \frac{1}{2} \mathbf{U}^T \left( \sum_{e=1}^{N_e} \mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e \right) \mathbf{U}
$$

The total strain energy can also be expressed directly in terms of the global stiffness matrix $\mathbf{K}$:

$$
U_{\text{total}} = \frac{1}{2} \mathbf{U}^T \mathbf{K} \mathbf{U}
$$

Comparing the two expressions yields the fundamental equation for global stiffness assembly:

$$
\mathbf{K} = \sum_{e=1}^{N_e} \mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e
$$

### Concrete example

Let us ground this with a highly specific example. Consider a 1D bar partitioned into 2 elements (3 nodes total: nodes 1, 2, 3), where each node has 1 DOF. The global system has $N=3$ DOFs.
Element 1 connects nodes 1 and 2. Element 2 connects nodes 2 and 3.

The Boolean matrix for Element 1, $\mathbf{L}_1$, maps the $2$ local DOFs of the element to the $3$ global DOFs:

$$
\mathbf{L}_1 = 
\begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0
\end{bmatrix}
$$

The transpose $\mathbf{L}_1^T$ (the "scatter" matrix) is:

$$
\mathbf{L}_1^T = 
\begin{bmatrix}
1 & 0 \\
0 & 1 \\
0 & 0
\end{bmatrix}
$$

Let the element stiffness matrix be:

$$
\mathbf{K}_1 = 
\begin{bmatrix}
k_{11} & k_{12} \\
k_{21} & k_{22}
\end{bmatrix}
$$

Now, evaluate the matrix product $\mathbf{L}_1^T \mathbf{K}_1 \mathbf{L}_1$:

1.  **First multiplication $\mathbf{K}_1 \mathbf{L}_1$:**

$$
\mathbf{K}_1 \mathbf{L}_1 = 
\begin{bmatrix}
k_{11} & k_{12} \\
k_{21} & k_{22}
\end{bmatrix}
\begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0
\end{bmatrix}
=
\begin{bmatrix}
k_{11} & k_{12} & 0 \\
k_{21} & k_{22} & 0
\end{bmatrix}
$$

2.  **Second multiplication $\mathbf{L}_1^T (\mathbf{K}_1 \mathbf{L}_1)$:**

$$
\mathbf{L}_1^T (\mathbf{K}_1 \mathbf{L}_1) = 
\begin{bmatrix}
1 & 0 \\
0 & 1 \\
0 & 0
\end{bmatrix}
\begin{bmatrix}
k_{11} & k_{12} & 0 \\
k_{21} & k_{22} & 0
\end{bmatrix}
=
\begin{bmatrix}
k_{11} & k_{12} & 0 \\
k_{21} & k_{22} & 0 \\
0 & 0 & 0
\end{bmatrix}
$$

The operation $\mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e$ mathematically expands a small $n_e \times n_e$ matrix into an $N \times N$ matrix, padding it with zeros everywhere except at the specific rows and columns corresponding to the element's global DOFs. The sum $\sum \mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e$ then superimposes all these padded matrices.

---

## Python implementation with `numpy.ix_`

While the equations above are mathematically elegant, constructing large, mostly empty $\mathbf{L}_e$ matrices and performing the full matrix multiplications is computationally intractable ($\mathcal{O}(N^3)$ operations for mostly zeros). In practice, the operation $\mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e$ is realized as an array slicing operation: "scatter the entries of $\mathbf{K}_e$ into the appropriate rows and columns of $\mathbf{K}$".

In Python, this is executed using `numpy.ix_`.

### Why `numpy.ix_` is the right tool

Advanced indexing in NumPy allows you to select arbitrary items based on their N-dimensional index. If you have an element whose local DOFs map to global DOFs `[0, 2]`, you want to add the $2 \times 2$ matrix $\mathbf{K}_e$ to rows `[0, 2]` and columns `[0, 2]` of $\mathbf{K}$.

If you simply write `K[[0, 2], [0, 2]]`, NumPy performs *integer array indexing*, pairing the indices up: it selects `K[0, 0]` and `K[2, 2]`, resulting in a 1D array of two elements. This is **not** a submatrix (block) selection!

To select the $2 \times 2$ block formed by the intersection of rows `[0, 2]` and columns `[0, 2]`, you need an open mesh grid. `numpy.ix_` takes 1D sequences and returns a tuple of N-dimensional arrays that broadcast to form the full block.

```python
import numpy as np

indices = [0, 2]
ix_grid = np.ix_(indices, indices)

print(ix_grid)
# Output:
# (array([[0],
#         [2]]), 
#  array([[0, 2]]))
```

Notice the shapes:
*   The row index array has shape `(2, 1)`: `[[0], [2]]`
*   The column index array has shape `(1, 2)`: `[[0, 2]]`

When NumPy uses these arrays to index `K`, the shapes broadcast together to form a `(2, 2)` grid:
*   Row 0 with Col 0 $\rightarrow (0, 0)$
*   Row 0 with Col 2 $\rightarrow (0, 2)$
*   Row 2 with Col 0 $\rightarrow (2, 0)$
*   Row 2 with Col 2 $\rightarrow (2, 2)$

This precisely targets the block entries corresponding to the selected global DOFs and matches the mathematical result of $\mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e$.

---

## Sparse and dense assembly in `assmk`

When accumulating $\mathbf{K}_e$ into $\mathbf{K}$, the syntax differs slightly based on whether $\mathbf{K}$ is a dense `numpy.ndarray` or a `scipy.sparse.lil_matrix`. 

### Dense arrays

A dense matrix stores every coefficient in continuous memory blocks. In Python, an in-place addition to a block using advanced indexing works natively:

```python
K[np.ix_(indices, indices)] += Ke
```

This invokes the C-level `__iadd__` routine. Because the memory is contiguous and fully allocated, NumPy computes the exact memory offsets and rapidly adds the values in-place without creating intermediate large array copies. It is blazingly fast for small meshes (e.g., $N < 1000$).

### SciPy sparse matrices

Finite element matrices are inherently sparse. For large systems ($N \gg 1000$), storing a dense matrix consumes $\mathcal{O}(N^2)$ memory, which quickly exceeds available RAM. Sparse matrices store only non-zero entries, requiring $\mathcal{O}(N)$ memory.

During assembly, SciPy's `lil_matrix` (List of Lists) is typically used. In this format, each row is a Python list of column indices and a corresponding list of values. 

For sparse arrays, augmented assignment with advanced indexing is problematic. The expression:

```python
K[np.ix_(indices, indices)] += Ke
```
fails or is highly inefficient for `lil_matrix` because `+=` attempts to mutate the extracted sparse slice in-place, which isn't robustly supported by SciPy's sparse advanced indexing mechanics. 

Instead, we must write:

```python
K[np.ix_(indices, indices)] = K[np.ix_(indices, indices)] + Ke
```

**What happens here step-by-step?**
1.  **Extract:** `K[np.ix_(indices, indices)]` extracts the sparse block, instantiating it temporarily as a dense or sparse sub-matrix.
2.  **Add:** `+ Ke` adds the dense element matrix to the extracted block.
3.  **Assign:** The `=` operator calls `__setitem__` on the `lil_matrix`. The LIL matrix updates its internal lists. If a zero entry suddenly becomes non-zero, it dynamically appends the new column index and value to the corresponding row list. 

While allocating lists during assignment has some overhead, it prevents the quadratic memory blow-up of a dense matrix and keeps the workflow practical for large meshes. After assembly is complete, the `lil_matrix` is typically converted to Compressed Sparse Row (`csr_matrix`) format for fast solving.

---

## Complete example: 3-element manual assembly

Below is a self-contained, 50-line runnable script demonstrating everything discussed: Boolean matrices via `ix_`, dense vs. sparse logic, and a 3-element system (4 nodes, 1 DOF per node).

```python
import numpy as np
import scipy.sparse as sp

def assemble_system(nn, elements, Ke, use_sparse=False):
    """
    nn: Total number of global nodes
    elements: List of node-pair tuples, e.g., [(0, 1), (1, 2), (2, 3)]
    Ke: Element stiffness matrix (assumed identical for all elements)
    """
    # 1. Initialize the global stiffness matrix
    if use_sparse:
        K = sp.lil_matrix((nn, nn), dtype=float)
    else:
        K = np.zeros((nn, nn), dtype=float)
        
    # 2. Iterate through elements
    for element_nodes in elements:
        indices = list(element_nodes)
        idx_grid = np.ix_(indices, indices)
        
        # 3. Apply Scatter-Add Logic
        if use_sparse:
            # Sparse requires extraction, addition, and reassignment
            K[idx_grid] = K[idx_grid] + Ke
        else:
            # Dense allows direct in-place addition
            K[idx_grid] += Ke
            
    return K

if __name__ == "__main__":
    # Define a 3-element bar connected in series: Node 0 -> 1 -> 2 -> 3
    nodes_total = 4
    topology = [(0, 1), (1, 2), (2, 3)]
    
    # Define a generic 2x2 local element stiffness matrix
    k_local = np.array([[ 100.0, -100.0],
                        [-100.0,  100.0]])
    
    # Execute Dense Assembly
    K_dense = assemble_system(nodes_total, topology, k_local, use_sparse=False)
    print("--- Dense Assembly Result ---")
    print(K_dense)
    print("\n--- Sparse Assembly Result (CSR format) ---")
    
    # Execute Sparse Assembly
    K_sparse_lil = assemble_system(nodes_total, topology, k_local, use_sparse=True)
    K_sparse_csr = K_sparse_lil.tocsr() # Convert to CSR for solving
    print(K_sparse_csr.toarray())
```
