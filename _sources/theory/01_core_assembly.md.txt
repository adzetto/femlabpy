# Global Matrix Assembly in Finite Element Analysis

## Introduction to Assembly

In the Finite Element Method (FEM), the behavior of a complex continuous domain is modeled by discretizing it into smaller, manageable subdomains called elements. Each element's behavior is characterized by local matrices—such as the element stiffness matrix $\mathbf{K}_e$ and the element internal force vector $\mathbf{q}_e$. To solve the system as a whole, these local contributions must be aggregated into a global system of equations representing the entire structure. This process is known as **global matrix assembly**.

The assembly procedure maps the local degrees of freedom (DOFs) of an individual element to the global degrees of freedom of the complete mesh structure. Mathematically, this mapping is elegantly represented using Boolean connectivity matrices. Computationally, however, using explicit Boolean matrix multiplication is highly inefficient. Instead, modern FEA codes, such as `femlabpy`, utilize direct index mapping and sophisticated memory management techniques (such as sparse matrix structures) to efficiently accumulate local contributions.

This chapter details both the rigorous mathematical foundation of assembly using Boolean matrices and the practical, high-performance implementation in Python using NumPy and SciPy.

---

## Mathematical Formulation: Boolean Assembly Matrices

The global stiffness matrix $\mathbf{K}$ and global internal force vector $\mathbf{q}$ are constructed by summing the contributions from all $N_e$ elements in the mesh. Let the global displacement vector be $\mathbf{U}$ and the local displacement vector for element $e$ be $\mathbf{u}_e$.

The relationship between the global degrees of freedom and the local degrees of freedom for an element $e$ is established through a Boolean connectivity matrix $\mathbf{L}_e$. The matrix $\mathbf{L}_e$ is a rectangular matrix of size $n_e \times N$, where $n_e$ is the number of local DOFs for the element and $N$ is the total number of global DOFs.

### Mapping Local to Global DOFs

The local displacement vector $\mathbf{u}_e$ is extracted from the global displacement vector $\mathbf{U}$ by:

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

Because each local DOF maps to exactly one global DOF, each row of $\mathbf{L}_e$ contains exactly one `1`, with all other entries being `0`. 

### Strain Energy and Global Stiffness

To derive the global assembly equation for the stiffness matrix, consider the total strain energy of the system, $U_{total}$, which is the sum of the strain energies of the individual elements:

$$
U_{total} = \sum_{e=1}^{N_e} U_e = \sum_{e=1}^{N_e} \frac{1}{2} \mathbf{u}_e^T \mathbf{K}_e \mathbf{u}_e
$$

Substitute the Boolean mapping $\mathbf{u}_e = \mathbf{L}_e \mathbf{U}$ into the energy equation:

$$
U_{total} = \sum_{e=1}^{N_e} \frac{1}{2} (\mathbf{L}_e \mathbf{U})^T \mathbf{K}_e (\mathbf{L}_e \mathbf{U}) = \frac{1}{2} \mathbf{U}^T \left( \sum_{e=1}^{N_e} \mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e \right) \mathbf{U}
$$

The total strain energy can also be expressed directly in terms of the global stiffness matrix $\mathbf{K}$:

$$
U_{total} = \frac{1}{2} \mathbf{U}^T \mathbf{K} \mathbf{U}
$$

Comparing the two expressions, we obtain the fundamental equation for global stiffness assembly:

$$
\mathbf{K} = \sum_{e=1}^{N_e} \mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e
$$

### Virtual Work and Global Forces

Similarly, the total internal virtual work $\delta W_{int}$ is the sum of the internal virtual work of each element. For a virtual displacement field $\delta \mathbf{U}$, the local virtual displacements are $\delta \mathbf{u}_e = \mathbf{L}_e \delta \mathbf{U}$.

$$
\delta W_{int} = \sum_{e=1}^{N_e} \delta \mathbf{u}_e^T \mathbf{q}_e = \sum_{e=1}^{N_e} (\mathbf{L}_e \delta \mathbf{U})^T \mathbf{q}_e = \delta \mathbf{U}^T \left( \sum_{e=1}^{N_e} \mathbf{L}_e^T \mathbf{q}_e \right)
$$

The global internal force vector $\mathbf{q}$ relates to the virtual work by $\delta W_{int} = \delta \mathbf{U}^T \mathbf{q}$. Thus, the assembly of the global internal force vector is:

$$
\mathbf{q} = \sum_{e=1}^{N_e} \mathbf{L}_e^T \mathbf{q}_e
$$

While the equations $\mathbf{K} = \sum \mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e$ and $\mathbf{q} = \sum \mathbf{L}_e^T \mathbf{q}_e$ are mathematically elegant, constructing large, mostly empty $\mathbf{L}_e$ matrices and performing the matrix multiplications is computationally intractable. In practice, the operation $\mathbf{L}_e^T \mathbf{K}_e \mathbf{L}_e$ simply means "scatter the entries of $\mathbf{K}_e$ into the appropriate rows and columns of $\mathbf{K}$". This scattering is implemented via direct array indexing.

---

## Python Implementation: The `assmk` and `assmq` Functions

In `femlabpy`, the theoretical Boolean mappings are translated into highly optimized array indexing operations using NumPy. The core functions handling this are `assmk` (assemble stiffness) and `assmq` (assemble forces), located in `src/femlabpy/assembly.py`.

### Extracting Global Degrees of Freedom

Before we can map $\mathbf{K}_e$ to $\mathbf{K}$, we must determine the global indices corresponding to the element's local DOFs. The topology vector `Te` dictates which nodes make up the element. For an element with nodes $[n_1, n_2, \dots]$, and a specified number of `dof` per node, the global indices are calculated. For example, if $n_i$ is 0-indexed and `dof = 2`, node $n_i$ maps to global DOFs $2n_i$ and $2n_i + 1$.

```python
from femlabpy._helpers import topology_nodes, node_dof_indices

# Extract the node IDs from the topology array (ignoring property IDs)
element_nodes = topology_nodes(Te)

# Compute the global DOF indices for these nodes
indices = node_dof_indices(element_nodes, dof)
```

### Assembling the Stiffness Matrix (`assmk`)

The `assmk` function takes the global stiffness matrix `K`, the element stiffness matrix `Ke`, the element topology `Te`, and the `dof` per node. It uses NumPy's `np.ix_` to cleanly define an open mesh grid of row and column indices. 

The `np.ix_` function constructs index arrays that broadcast properly to access the block in the matrix corresponding to `indices` for both rows and columns.

```python
def assmk(K, Ke, Te, dof: int = 1):
    element_nodes = topology_nodes(Te)
    indices = node_dof_indices(element_nodes, dof)
    element_matrix = as_float_array(Ke)
    
    if is_sparse(K):
        # SciPy sparse matrices (like LIL) do not support augmented assignment `+=` 
        # seamlessly with advanced indexing in the same way dense arrays do.
        K[np.ix_(indices, indices)] = K[np.ix_(indices, indices)] + element_matrix
    else:
        # For dense NumPy arrays, in-place addition is highly optimized.
        K[np.ix_(indices, indices)] += element_matrix
    return K
```

Notice the conditional check `is_sparse(K)`. Advanced slicing works differently between standard dense `numpy.ndarray` and `scipy.sparse` matrix formats, necessitating the split logic.

### Assembling the Internal Force Vector (`assmq`)

The `assmq` function performs the analog operation for a 1D array. Since vectors only require one-dimensional indexing, `np.ix_` is unnecessary. We directly index the rows of `q`.

```python
def assmq(q, qe, Te, dof: int = 1):
    element_nodes = topology_nodes(Te)
    indices = node_dof_indices(element_nodes, dof)
    
    q = as_float_array(q)
    qe = as_float_array(qe).reshape(-1, 1)
    
    # In-place accumulation of internal forces
    q[indices, 0] += qe[:, 0]
    return q
```

The reshapes guarantee that `qe` acts mathematically like a column vector, maintaining consistent behavior across operations.

---

## Performance Considerations: Dense vs. Sparse Assembly

When initializing the system via `femlabpy.core.init`, users can specify whether to use dense arrays or sparse matrices. 

```python
def init(nn: int, dof: int, *, dynamic: bool = False, use_sparse: bool | None = None):
    total_dofs = int(nn) * int(dof)
    if use_sparse is None:
        use_sparse = nn >= 1000  # Default heuristic
    # ...
```

By default, the code heuristically switches to sparse matrices when the number of nodes `nn >= 1000`. This distinction is vital due to computational complexity.

### Dense Arrays

A dense global stiffness matrix stores every coefficient, requiring $\mathcal{O}(N^2)$ memory, where $N$ is the total number of DOFs. For a problem with 10,000 DOFs, a dense matrix requires approximately $800$ MB of RAM ($10,000 \times 10,000 \times 8 \text{ bytes}$). For 100,000 DOFs, memory balloons to 80 GB, which is prohibitive for standard workstations.

Dense assembly using `K[np.ix_(indices, indices)] += element_matrix` is blazingly fast in Python because it routes directly to compiled C memory routines without overhead. For small systems (e.g., academic problems with fewer than 1000 nodes), the speed of pure dense array updates outweighs the $\mathcal{O}(N^2)$ memory cost. 

### Sparse Matrices (LIL / CSR format)

A finite element stiffness matrix is inherently sparse because elements only connect to adjacent elements. Nodes that do not share an element have a mutual stiffness coefficient of precisely zero. For typical finite element meshes, the number of non-zero entries per row is bounded by a small constant $C \ll N$ (usually $C < 100$), regardless of the mesh size. Thus, a sparse matrix requires only $\mathcal{O}(N)$ memory.

In `femlabpy`, if `use_sparse=True`, `zeros_matrix` generally initializes a `scipy.sparse.lil_matrix` (List of Lists format). 

1. **LIL Format during Assembly:** The LIL format is designed for efficient modification of sparsity structure (inserting non-zero elements). Under the hood, LIL matrices store row data as Python lists. When `assmk` executes:
   ```python
   K[np.ix_(indices, indices)] = K[np.ix_(indices, indices)] + element_matrix
   ```
   the matrix rapidly dynamically allocates memory for new non-zero elements.
   
2. **CSR Format during Solving:** Once assembly is complete and the matrix equation $\mathbf{K} \mathbf{U} = \mathbf{F}$ is ready to be solved, LIL matrices are typically converted to Compressed Sparse Row (CSR) or Compressed Sparse Column (CSC) formats. CSR formats are highly optimized for matrix-vector multiplications and sparse linear solvers (e.g., SuperLU), but are inefficient for assembly since inserting a new entry requires shifting large portions of contiguous memory.

The `assmk` logic safely abstracts away this complexity, simply providing the mathematical interface required to accumulate elemental stiffness contributions properly.

---

## Complete Example

Below is a runnable Python snippet demonstrating the theoretical assembly process applied to a small two-element truss system in 1D.

```python
import numpy as np
from femlabpy.core import init
from femlabpy.assembly import assmk, assmq

# 1D Truss: 3 Nodes, 2 Elements, 1 DOF per node
# Node 0 --- (El 1) --- Node 1 --- (El 2) --- Node 2

# Element Stiffness Matrices (EA/L = 100)
Ke1 = np.array([[ 100, -100],
                [-100,  100]])
Ke2 = np.array([[ 100, -100],
                [-100,  100]])

# Internal Force Vectors (uniform load example)
qe1 = np.array([[5], [5]])
qe2 = np.array([[5], [5]])

# Topology definitions (Nodes, Property ID)
Te1 = np.array([0, 1, 99]) # Nodes 0, 1
Te2 = np.array([1, 2, 99]) # Nodes 1, 2

# 1. Initialize Global Matrices
K, p, q = init(nn=3, dof=1, use_sparse=False)

# 2. Assemble Element 1
K = assmk(K, Ke1, Te1, dof=1)
q = assmq(q, qe1, Te1, dof=1)

# 3. Assemble Element 2
K = assmk(K, Ke2, Te2, dof=1)
q = assmq(q, qe2, Te2, dof=1)

print("Global Stiffness Matrix K:")
print(K)
print("\nGlobal Internal Force Vector q:")
print(q)
```

**Output:**
```
Global Stiffness Matrix K:
[[ 100. -100.    0.]
 [-100.  200. -100.]
 [   0. -100.  100.]]

Global Internal Force Vector q:
[[ 5.]
 [10.]
 [ 5.]]
```

In the output, note how global DOF `1` (the middle node) naturally accumulates stiffness ($100 + 100 = 200$) and forces ($5 + 5 = 10$) from both adjacent elements, perfectly capturing the physics of connectivity through the logic of Boolean matrices realized via NumPy indexing.
