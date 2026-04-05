# Chapter 8: Custom Element Development

`femlabpy` is built around small, explicit kernels instead of a deep class hierarchy. That makes custom element work practical: you write the element math once, then wire it into the same assembly and recovery flow used everywhere else in the library.

The important part is not just the local matrix formula. A usable custom element has to match the repository conventions for topology rows, material tables, degrees of freedom, and recovery outputs. This chapter shows that flow from the helper functions up to the assembled global system.

## 8.1 The element contract

In this repository, every element family follows the same pattern:

`ke...`
: Return a local stiffness matrix.

`q...`
: Return the local internal-force vector and, when useful, element output quantities such as stress and strain.

`m...`
: Return a local mass matrix for dynamic analysis.

`k...`
: Assemble all element stiffness matrices into the global matrix.

`q...`
: Assemble all element internal-force vectors into the global vector.

That split is deliberate. The local routine should know only element geometry, materials, and interpolation. The assembly routine should know only how to loop through the topology table and scatter the element result into the global system.

The helper functions in `src/femlabpy/_helpers.py` define the data model that makes this work:

`topology_nodes`
: Removes the last entry of a legacy topology row and returns the connected node numbers.

`topology_property`
: Reads the final entry of the same topology row and treats it as the material or property id.

`material_row`
: Converts that one-based property id into the matching row from the material table.

`node_dof_indices`
: Expands one-based node numbers into flattened global DOF indices.

If you keep those four conventions intact, your custom element can plug into the rest of the package without special-case glue.

## 8.2 How assembly works here

The assembly helpers in `src/femlabpy/assembly.py` are intentionally small:

`assmk(K, Ke, Te, dof)`
: Extract the node numbers from `Te`, build the global DOF indices, and add `Ke` into the global matrix `K`.

`assmq(q, qe, Te, dof)`
: Do the same thing for a force vector.

That means the assembly loop in a custom driver should stay boring. Boring is good here. The real work belongs in the local kernel.

For a dense matrix, the pattern is straightforward:

```python
indices = node_dof_indices(topology_nodes(row), dof)
K[np.ix_(indices, indices)] += Ke
```

For a sparse matrix, the same indexing path is used, but you should keep the matrix in an assignment-friendly sparse format during assembly and convert it later if needed.

The practical implication is simple:

1. Build local matrices with NumPy.
2. Assemble with `assmk` and `assmq`.
3. Keep the topology format consistent with the rest of the library.
4. Do not mix local kernel logic with global solve logic.

## 8.3 Scalar potential elements

The scalar potential family is the cleanest place to start because it uses `dof = 1`. The repository already includes the full pattern in `src/femlabpy/elements/triangles.py` and `src/femlabpy/elements/quads.py`.

The triangular scalar kernel is `ket3p`:

```python
def ket3p(Xe, Ge):
    a, area = _triangle_geometry(Xe)
    props = as_float_array(Ge).reshape(-1)
    conductivity = props[0]
    D = np.eye(2, dtype=float) * conductivity
    B = (1.0 / (2.0 * area)) * np.column_stack([-a[:, 1], a[:, 0]]).T
    Ke = area * B.T @ D @ B
    if props.size > 1:
        b = props[1]
        Ke = Ke + (b * area / 12.0) * np.array(
            [[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]]
        )
    return Ke
```

What matters in that routine:

`Xe`
: The element coordinates. The kernel never asks where the element came from; it only needs the coordinates of the current triangle.

`Ge`
: A material row. For this element, the first value is conductivity and the optional second value is a reaction term.

`B`
: The gradient operator. It maps nodal potentials to spatial gradients.

`Ke`
: The conductivity matrix that will be scattered into the global system.

The companion recovery kernel is `qet3p`, which returns the equivalent nodal flux vector plus the element gradient and flux:

```python
qe, Se, Ee = qet3p(Xe, Ge, Ue)
```

`qt3p` then loops over the topology table, calls `qet3p` for each element, and accumulates the global vector with `assmq`.

The quadrilateral scalar path follows the same structure. `kq4p` uses the same topological and assembly conventions, but the local integration rule is a 2x2 Gauss scheme instead of the triangle's closed-form expression. That difference stays inside the element kernel, not the public API. If you understand one family, the other family is the same workflow with different interpolation.

## 8.4 Structural elements

The structural triangle kernels use the same assembly model, but the local DOF layout changes to `dof = 2`.

`ket3e`
: Returns the 6x6 plane-stress or plane-strain stiffness matrix for a constant-strain triangle.

`qet3e`
: Returns the internal-force vector together with element stress and strain values.

`kt3e`
: Assembles all T3 element stiffness matrices into the global stiffness matrix.

`qt3e`
: Recomputes element stress and strain from the global displacement field and scatters the equivalent nodal forces.

The ordering is:

```text
[u1, v1, u2, v2, u3, v3]
```

That ordering is what `element_dof_indices` and `assmk` expect. If your custom element uses a different ordering, every downstream routine will be wrong even if the local math is correct.

The same pattern appears in the 3D solid family:

`keT4e`, `qeT4e`, `kT4e`, `qT4e`
: Four-node tetrahedron element kernels and assembly wrappers.

`keh8e`, `qeh8e`, `kh8e`, `qh8e`
: Eight-node hexahedron element kernels and assembly wrappers.

These modules are useful references because they show how the repository handles:

1. Vectorized batch assembly.
2. Material selection through topology property ids.
3. Stress and strain recovery stored per element.
4. Sparse-aware accumulation without changing the public API.

## 8.5 Writing a new kernel

When you add a new element, keep the implementation split in the same way as the existing code.

1. Write a local stiffness function `kefoo(Xe, Ge)` that takes only one element's coordinates and material data.
2. Write a recovery function `qfoo(Xe, Ge, Ue)` if the element has internal forces, stresses, or fluxes.
3. Write a global assembly wrapper `kfoo(K, T, X, G)` that loops over `T`, extracts `Xe` with `topology_nodes`, resolves the material row with `material_row`, and calls `assmk`.
4. Write the matching force wrapper `qfoo_global(q, T, X, G, u)` and call `assmq`.
5. Add the symbols to `src/femlabpy/elements/__init__.py` and, if needed, to the package root exports.

The following skeleton matches the style used in the repository:

```python
def kefoo(Xe, Ge):
    Xe = as_float_array(Xe)
    props = as_float_array(Ge).reshape(-1)
    # build local interpolation, gradients, constitutive matrix
    # return the element matrix only
    return Ke


def qfoo(Xe, Ge, Ue):
    Xe = as_float_array(Xe)
    Ue = as_float_array(Ue).reshape(-1, 1)
    # compute local response
    return qe, Se, Ee


def kfoo(K, T, X, G):
    coords = as_float_array(X)
    dof = 2  # change to 1 for scalar kernels
    for row in as_float_array(T):
        nodes = topology_nodes(row)
        prop = topology_property(row)
        Ke = kefoo(coords[nodes - 1], material_row(G, prop))
        K = assmk(K, Ke, row, dof)
    return K


def qfoo_global(q, T, X, G, u):
    coords = as_float_array(X)
    dof = 2  # change to 1 for scalar kernels
    U = as_float_array(u).reshape(coords.shape[0], dof)
    for row in as_float_array(T):
        nodes = topology_nodes(row)
        prop = topology_property(row)
        Ue = U[nodes - 1].reshape(-1, 1)
        qe, _, _ = qfoo(coords[nodes - 1], material_row(G, prop), Ue)
        q = assmq(q, qe, row, dof)
    return q
```

Do not try to make the local kernel aware of the global matrix shape. That coupling makes the code harder to test and much harder to reuse.

## 8.6 Example integration checklist

Before a new element is considered finished, check these items:

1. The local kernel returns the expected matrix shape.
2. The recovery kernel returns arrays with predictable shapes and consistent ordering.
3. The topology row format matches the rest of the package.
4. The material table row format is documented and used consistently.
5. Dense and sparse assembly both work on a small test problem.
6. The element is exported from `src/femlabpy/elements/__init__.py`.
7. The element is reachable from the package root if users are meant to import it directly.
8. A short documentation example shows the full flow from `Xe` and `Ge` to assembled `K` and `q`.

A regression test should use a tiny problem where you can compute the same result by hand. For a new scalar element, a single triangle or quadrilateral is enough. For a structural element, compare the new kernel against a known textbook case or one of the existing `femlabpy` elements with the same interpolation order.

## 8.7 What to read next

If you are adding a new kernel family, read these files in this order:

`src/femlabpy/_helpers.py`
: For the node, material, and index conventions.

`src/femlabpy/assembly.py`
: For the scatter-add pattern.

`src/femlabpy/elements/triangles.py`
: For the simplest full examples in both scalar and structural form.

`src/femlabpy/elements/quads.py`
: For the higher-order version of the same workflow.

`src/femlabpy/elements/solids.py`
: For vectorized batch assembly and 3D indexing.

If you keep the math local and the assembly generic, a custom element stays small enough to understand and test.
