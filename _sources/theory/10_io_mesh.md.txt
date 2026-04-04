# Mesh I/O, Data Types, and Helpers

In computational mechanics, finite element analysis (FEA) heavily relies on robust and efficient representation of mesh data and sparse linear algebra operations. The `femlabpy` framework is designed to bridge the gap between traditional, legacy MATLAB implementations and modern, Pythonic object-oriented paradigms. 

This chapter provides a comprehensive overview of how mesh topologies are imported, the data structures used to store them, and the underlying mathematical and programmatic helpers that facilitate the assembly and solution of the resulting finite element systems.

---

## 1. The Gmsh File Format (`.msh`)

The open-source mesh generator **Gmsh** is the standard source of finite element meshes for `femlabpy`. Over time, Gmsh has evolved its `.msh` file format. `femlabpy` must accommodate both the legacy 2.2 format and the modern 4.1 format.

### 1.1 The Legacy v2.2 Format

The Gmsh v2.2 format is characterized by its simple ASCII block structure, which is easily parsed by line-by-line readers. The file is divided into distinct sections such as `$MeshFormat`, `$PhysicalNames`, `$Nodes`, and `$Elements`.

A typical v2.2 node section looks like this:
```text
$Nodes
number_of_nodes
node_number x-coord y-coord z-coord
...
$EndNodes
```
The coordinate vector for a node $i$ is denoted as $\mathbf{x}_i \in \mathbb{R}^3$.

The element section defines the topological connectivity:
```text
$Elements
number_of_elements
element_number element_type number_of_tags <tags> node_number_list
...
$EndElements
```
Here, the `element_type` corresponds to standardized integer IDs (e.g., `1` for lines, `2` for triangles, `3` for quadrilaterals, `4` for tetrahedrons). The `tags` usually encapsulate the physical and elementary regions, representing the material property IDs or boundary condition labels.

### 1.2 The Modern v4.1 Format

The modern Gmsh v4.1 format introduces block-based data storage, which significantly improves read/write performance for extremely large meshes. However, the block-based layout is incompatible with simple line-by-line procedural parsers written for v2.2.

To natively support v4.1 without rewriting the entire legacy logic, `femlabpy` employs a strategic workaround: if the Python `gmsh` SDK is installed, `femlabpy` loads the v4.1 mesh into the SDK's internal memory and temporarily re-emits it as a v2.2 legacy ASCII file. This guarantees backward compatibility and unified parsing behavior:

```python
import gmsh
gmsh.initialize(readConfigFiles=False)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.open("mesh_v41.msh")
gmsh.write("legacy_v22_ascii.msh")
```

---

## 2. Parsing Meshes with `load_gmsh2`

The `load_gmsh2` parser in `femlabpy.io.gmsh` is responsible for translating the `.msh` text structure into structured numpy arrays.

### 2.1 Physical Groups and Topology Arrays

A critical step in FEA preprocessing is associating specific regions of the mesh with physical properties (e.g., Young's modulus, mass density). In Gmsh, these are defined via "Physical Groups".

When `load_gmsh2` processes an element row, it extracts the tags. In the v2.2 format, the first tag in the `<tags>` list corresponds to the Physical Group ID. 

```python
def load_gmsh2(filename, which=None) -> GmshMesh:
    # ...
    # Within the element parsing loop:
    # tags = element_parts[3 : 3 + num_tags]
    # The first tag is typically the Physical Group ID
```

This Physical Group ID directly maps to the material property IDs in the topology array `T`. For example, in a simple 2D triangle element, `T` is a row matrix where the first 3 columns are the node numbers, and the last column is the material property ID:
$$ T_e = [n_1, n_2, n_3, \text{PhysicalGroupID}] $$

`femlabpy` organizes elements by type and appends this tag to the explicit topology arrays, enabling quick lookup of material properties during integration.

---

## 3. The `GmshMesh` Data Structure

To manage the parsed data, `femlabpy.types` implements the `GmshMesh` dataclass. This structure bridges the gap between modern Python attributes and the legacy MATLAB variables.

### 3.1 Nodal and Element Arrays

At its core, `GmshMesh` maintains the nodal positions $\mathbf{X} \in \mathbb{R}^{N_{nodes} \times 3}$, and aggregated arrays like `element_infos` and `element_tags`.

```python
@dataclass
class GmshMesh:
    positions: np.ndarray          # Shape: (N, 3)
    element_infos: np.ndarray      # Element ID, Type, Tag count, Dimension
    element_tags: np.ndarray       # Tag values
    element_nodes: np.ndarray      # Connectivity
    # ...
```

### 3.2 Explicit Element Topologies

The class provides explicit, zero-padded numpy arrays for different element types:
- `points`: 1-node points (Type 15)
- `lines`: 2-node lines (Type 1)
- `triangles`: 3-node triangles (Type 2)
- `quads`: 4-node quadrilaterals (Type 3)
- `tets`: 4-node tetrahedrons (Type 4)
- `hexa`: 8-node hexahedrons (Type 5)

Each of these arrays holds the node numbers followed by the physical tag. For instance, the `triangles` array has a shape of `(N_triangles, 4)`.

### 3.3 Legacy Aliases

To maintain backward compatibility with classroom scripts that expect MATLAB variable names, `GmshMesh` uses a custom `__getattr__` implementation to provide upper-case aliases:
- `mesh.POS` resolves to `mesh.positions`
- `mesh.TRIANGLES` resolves to `mesh.triangles`
- `mesh.nbTriangles` resolves to the number of rows in `mesh.triangles`

---

## 4. Helper Functions for Computational Mechanics

Once the mesh is loaded, `femlabpy._helpers` provides the essential building blocks for matrix assembly, degree-of-freedom mapping, and linear system solving.

### 4.1 Array Coercion

In Python, sequences can take many forms (lists, tuples, raw numpy arrays). FEA routines require strict dimensional and datatype constraints.

**`as_float_array`**: Coerces any array-like structure into a floating-point numpy array. It uses `dtype=float` internally.
```python
from femlabpy._helpers import as_float_array
node_coords = as_float_array([[0, 0], [1, 0], [0, 1]])
```

**`as_column`**: Specifically reshapes the input to a 2D column vector $\mathbf{v} \in \mathbb{R}^{n \times 1}$, which is critical when assembling force vectors $\mathbf{F}$.

### 4.2 Matrix Operations and Sparsity

FEA generates large, highly sparse stiffness matrices $\mathbf{K}$. Knowing whether a matrix is dense or sparse dictates how operations like matrix factorization and solution are handled.

**`is_sparse`**: Uses `scipy.sparse.issparse` to verify if a given matrix is sparse.
```python
from femlabpy._helpers import is_sparse, zeros_matrix
K = zeros_matrix(100, use_sparse=True)
print(is_sparse(K))  # Returns: True
```

### 4.3 Solving Linear Systems

The core of any finite element code is solving the algebraic system:
$$ \mathbf{K} \mathbf{U} = \mathbf{F} $$
where $\mathbf{K}$ is the global stiffness matrix, $\mathbf{U}$ is the vector of nodal displacements, and $\mathbf{F}$ is the global force vector.

`femlabpy` provides the **`solve_linear_system`** function to handle this generically. It checks the sparsity of $\mathbf{K}$ and routes the calculation to the most appropriate solver:
- For **sparse** matrices, it uses `scipy.sparse.linalg.spsolve`.
- For **dense** matrices, it uses `numpy.linalg.solve`.

```python
import numpy as np
from femlabpy._helpers import solve_linear_system, zeros_vector
import scipy.sparse as sp

# Assemble a simple system
K_sparse = sp.eye(3, format="csc") * 2.0
F = as_column([10.0, 20.0, 30.0])

# Solve the system
U = solve_linear_system(K_sparse, F)
print(U.flatten()) # Outputs: [ 5. 10. 15.]
```

For tiny classroom examples, the legacy MATLAB behavior involved fallback heuristics based on the condition number of the system. `femlabpy` implements `solve_legacy_symmetric_system` which gracefully falls back to `scipy.linalg.solve` with `assume_a='pos'` (Cholesky decomposition) or `'sym'` (symmetric indefinite factorization) depending on the evaluated matrix condition.

---

## 5. Summary

The infrastructure laid out in `femlabpy` perfectly encapsulates the raw strings of a `.msh` file into highly optimized, algebraically-ready array structures. By leveraging `GmshMesh` and the robust numeric helpers, engineers can abstract away file I/O and topological parsing, focusing solely on the assembly of the stiffness operators and the mathematical formulations of computational mechanics.
