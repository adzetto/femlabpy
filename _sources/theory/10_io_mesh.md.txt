# Mesh I/O and Gmsh Import

This chapter explains how `femlabpy` mathematically and structurally turns a `.msh` file into the primary array structures used by the element kernels: the coordinate matrix $\mathbf{X}$ and the topology matrix $\mathbf{T}$.

## 1. The Core Data Structures

In `femlabpy`, the finite element geometry is entirely defined by two globally ordered matrices:

1.  **Coordinate Matrix ($\mathbf{X}$):** A floating-point array of shape `(nn, ndim)` where `nn` is the total number of nodes and `ndim` is the spatial dimension (2 or 3). The $i$-th row $\mathbf{X}[i, :]$ contains the $[x, y, z]$ coordinates of the $i$-th node.
2.  **Topology Matrix ($\mathbf{T}$):** An integer array representing element connectivity. For a specific element type with $k$ nodes (e.g., $k=4$ for Q4 elements), $\mathbf{T}$ has shape `(ne, k+1)` where `ne` is the number of elements of that type. The first $k$ columns contain the zero-based node indices that make up the element, and the final column contains the physical property tag/material ID assigned to that element.

## 2. Parsing the Mesh File

`femlabpy` supports both Gmsh 2.2 ASCII layout and modern 4.x binary/ASCII layouts. Modern meshes are converted into the 2.2 legacy layout in memory using the `gmsh` Python SDK if available.

### 2.1. Nodal Transformation ($\$Nodes \to \mathbf{X}$)

In the Gmsh file, the `$Nodes` block defines the mesh nodes:

```text
$Nodes
nn
node_id_1  x_1  y_1  z_1
node_id_2  x_2  y_2  z_2
...
$EndNodes
```

**Transformation:**
Gmsh `node_id`s are generally 1-based and not necessarily contiguous. `femlabpy` reads these into a dense NumPy array $\mathbf{X}$. During this process, node IDs are normalized to a strict `0` to `nn-1` zero-based index system.

$$
\text{Gmsh Node } i \implies \mathbf{X}[i-1, :] = [x_i, y_i, z_i]
$$

### 2.2. Element Transformation ($\$Elements \to \mathbf{T}$)

The `$Elements` block defines the mesh connectivity and physical grouping:

```text
$Elements
ne
elm_id_1  elm_type  num_tags  tag_1  tag_2  ...  node_1  node_2  ...  node_k
...
$EndElements
```

**Transformation:**
1.  **Filtering:** Elements are filtered and grouped by their `elm_type` (e.g., type 2 = 3-node triangle, type 3 = 4-node quad).
2.  **Node Index Normalization:** The 1-based `node_i` values are decremented by 1 to match the normalized row indices of $\mathbf{X}$.
3.  **Property Tagging:** Gmsh allows multiple tags per element. `femlabpy` extracts the *first* physical tag (`tag_1`) and treats it as the material/property ID index.
4.  **Assembly:** For a specific element group, the topology matrix $\mathbf{T}$ is assembled such that row $e$ is:

$$
\mathbf{T}[e, :] = [ (\text{node}_1 - 1), (\text{node}_2 - 1), \dots, (\text{node}_k - 1), \text{tag}_1 ]
$$

This structure ensures that when iterating over elements, `femlabpy` can simultaneously extract the geometry via `Xe = X[T[e, :-1], :]` and the material row via `Ge = G[T[e, -1] - 1, :]`.

## 3. Loader Functions

`femlabpy.io.gmsh` exposes two public loaders that construct the `GmshMesh` container object containing these $\mathbf{X}$ and $\mathbf{T}$ matrices:

### 3.1. `load_gmsh`
`load_gmsh(filename)` reproduces the legacy `load_gmsh.m` semantics. It loads all explicit element tables for all supported element families.

### 3.2. `load_gmsh2`
`load_gmsh2(filename, which=None)` allows selective loading. You can pass a list of element types `which=[2, 3]` to only extract topology matrices for triangles and quads, saving memory for large meshes containing unused bounding elements.

## 4. Practical Extraction

Once loaded, the `GmshMesh` object exposes the parsed matrices natively:

```python
import femlabpy as fp

mesh = fp.load_gmsh2("model.msh", which=[3]) # Load only Q4 elements

# The coordinate matrix X (nn x 3)
X = mesh.positions

# The topology matrix T for quads (ne x 5)
T_q4 = mesh.quads 
```

This clean separation ensures that `femlabpy`'s assembly routines (like `assmk()`) remain completely decoupled from the file format, operating purely on dense, vectorized NumPy arrays.