# Mesh I/O and Gmsh Import

This chapter explains how `femlabpy` turns a `.msh` file into the array
structures used by the rest of the library. The important part is not the file
format itself. The important part is how node coordinates, element topology, and
physical tags become NumPy arrays that the element kernels can consume directly.

## Mesh Formats

`femlabpy` supports the older Gmsh 2.2 ASCII layout and the newer 4.x family.
The implementation prefers a single normalized parsing path instead of two
different loaders.

### Legacy ASCII meshes

Legacy 2.2 meshes are simple to parse because nodes and elements are stored in
plain text blocks. The parser reads the `$Nodes` block into a coordinate array
with shape `(nn, 3)` and the `$Elements` block into normalized element records.

The key convention is that Gmsh writes one-based node ids, while the Python
arrays inside `femlabpy` are zero-based. The loader resolves that mismatch as
part of parsing, so later assembly code can use direct NumPy indexing.

### Modern Gmsh meshes

Modern 4.x meshes use block-based storage. `femlabpy` handles those files by
using the optional official `gmsh` SDK when it is installed, then re-emitting
the mesh into the legacy ASCII layout expected by the parser. That keeps the
rest of the code path identical.

The practical effect is simple:

1. if the file is already compatible, it is parsed directly;
2. if the file is a modern mesh and the SDK is available, it is converted to a
   legacy ASCII view first;
3. if the SDK is missing, the loader fails early instead of guessing.

## Loader Functions

`femlabpy.io.gmsh` exposes two public loaders.

### `load_gmsh`

`load_gmsh(filename)` reproduces the legacy `load_gmsh.m` semantics. It loads
all explicit element tables for the supported element families and returns a
fully populated `GmshMesh` object.

Use this when you want the closest behavior to the original MATLAB loader and do
not mind carrying all explicit arrays.

### `load_gmsh2`

`load_gmsh2(filename, which=None)` is the more flexible loader. It returns the
same `GmshMesh` container, but it lets you choose which explicit element tables
should be materialized.

- `which=None` loads all explicit arrays.
- `which=-1` or an empty iterable skips the explicit arrays.
- a list such as `[2, 3, 4]` loads only the requested element types.

This matters when you only need one topology table, because the mesh object can
stay much smaller.

### What the loader returns

Both loaders return a `GmshMesh` object. The object keeps the normalized Python
fields used internally by `femlabpy` and the legacy aliases used by the old
FemLab scripts.

The most important fields are:

- `positions`: node coordinates, always stored as `(nn, 3)`.
- `element_infos`: a compact summary of each parsed element.
- `element_tags`: the element tag table.
- `element_nodes`: the element connectivity table.
- `nb_type`: counts per Gmsh element type.
- `bounds_min` and `bounds_max`: axis-aligned bounds of the mesh.

The explicit topology arrays are the ones most users care about:

- `triangles`, `quads`, `tets`, `hexa`, and the higher-order variants.
- each row stores the node ids for one element and the first element tag in the
  last column.

## GmshMesh Structure

`GmshMesh` is the normalization layer between the mesh file and the FEM code.
It is defined in `src/femlabpy/types.py` and is intentionally opinionated.

### Python fields

The lowercase fields are the stable API:

- `positions`
- `element_infos`
- `element_tags`
- `element_nodes`
- `nb_type`
- the explicit topology arrays such as `triangles` and `quads`

The class also stores `legacy_element_infos`, `legacy_element_tags`,
`loader_name`, `explicit_types`, and `nodes_per_type_of_element`. Those fields
help preserve the behavior of the original classroom loaders.

### Legacy aliases

The class implements `__getattr__` so older scripts can still use names such as
`POS`, `TRIANGLES`, `QUADS`, `nbTriangles`, `MIN`, and `MAX`.

That is not cosmetic. It lets the wrapper layer and the tutorial code keep the
old FemLab naming while the internals stay consistent with Python conventions.

## Parsing Flow

The private helpers in `io/gmsh.py` are worth understanding even if you never
call them directly.

### File version detection

`_mesh_format_version()` reads the `$MeshFormat` header and returns the declared
version when possible.

### Legacy view creation

`_legacy_view_path()` either returns the original file path or creates a
temporary 2.2-style ASCII view through the `gmsh` SDK. That is the compatibility
bridge for modern meshes.

### Normalized parsing

`_parse_gmsh_file()` converts the mesh into node arrays, bounds, and parsed
element dictionaries. `_build_normalized_mesh()` then turns those parsed records
into the final `GmshMesh` object.

The important point is that the parser keeps the connectivity and
the physical tags separate. That allows the explicit topology arrays to carry the
first physical tag in the last column, while the more general tables still keep
the full metadata.

## Practical Usage

```python
import femlabpy as fp

mesh = fp.load_gmsh2("model.msh", which=[2, 3, 4])

print(mesh.nbNod)
print(mesh.nbElm)
print(mesh.triangles.shape)
print(mesh.quads.shape)
print(mesh.bounds_min, mesh.bounds_max)
```

For code that needs to inspect the physical group associated with an element
row, use `mesh.property_numbers(...)` or read `mesh.element_tags` directly.

## Reading Checklist

If you are editing the loader, keep these rules in mind:

1. preserve one-based ids at the file boundary;
2. convert to zero-based NumPy indexing as soon as the data enters the mesh
   object;
3. keep the physical tag in the last column of the explicit topology arrays;
4. do not silently reinterpret malformed meshes.
