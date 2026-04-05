# Chapter 2: Element Library

The element library in `femlabpy` turns the continuous mechanics problem into local element matrices and vectors that can be assembled into a global system. The code follows a consistent pattern across all families:

- `Xe` holds one element's coordinates.
- `Ge` holds the compact material row used by that kernel.
- `T` holds the topology table for all elements.
- `K`, `M`, and `q` are the global stiffness, mass, and internal-force arrays.
- Functions starting with `k` return stiffness matrices, `m` return mass matrices, and `q` return internal forces, fluxes, stresses, or strains depending on the physics.

The most useful thing to keep in mind while reading the code is that the comments and docstrings are not decorative. They describe the exact array flow used in the implementation: how nodal coordinates are reshaped, how Jacobians are built, how `B` matrices are formed, and how each element contribution is scattered into the global matrix.

## 2.1 Common data layout

The library uses a compact data layout that mirrors the original FemLab conventions but is written in NumPy style.

- Element coordinates are always given as `Xe`, usually with shape `(nnode, ndim)`.
- Material rows are compact and element-specific. For example, a bar uses `[A, E, rho]`, a triangle may use `[E, nu]` or `[E, nu, type]`, and some quadrilateral kernels accept extra columns for thickness or history data.
- Topology rows are 1-based, because the code keeps the legacy MATLAB/FemLab style at the interface and subtracts 1 only when indexing NumPy arrays.
- Global vector ordering is nodal and block-based. In 2D vector problems the layout is `[u1, v1, u2, v2, ...]`. In scalar problems it is `[phi1, phi2, ...]`.

The kernel comments repeatedly use the same symbols:

- `J` or `Jt` is the Jacobian matrix.
- `B` is the strain-displacement matrix for mechanics or the gradient matrix for scalar problems.
- `D` or `C` is the constitutive matrix.
- `qe`, `qeT4e`, `qeq4e`, and similar names are element internal-force or flux vectors.
- `Se` and `Ee` usually store stress-like and strain-like quantities at integration points.

This naming is deliberate. It makes the code easy to scan once you learn the pattern.

## 2.2 Bar and truss elements

The bar family is the simplest place to see the library's design. The low-level routines are:

- `kebar(Xe0, Xe1, Ge)` for the tangent stiffness of a geometrically nonlinear bar.
- `qebar(Xe0, Xe1, Ge)` for the corresponding internal-force vector.
- `kbar(K, T, X, G, u=None)` for global stiffness assembly.
- `qbar(q, T, X, G, u=None)` for global internal-force assembly.
- `mebar(Xe, Ge, dof=2, lumped=False)` for one element mass matrix.
- `mbar(M, T, X, G, dof=2, lumped=False)` for global mass assembly.

The nonlinear bar kernel works with both the initial and current element geometry. That is why `kebar` and `qebar` take `Xe0` and `Xe1`. The code first computes the initial axis vector `a0`, the current axis vector `a1`, and the corresponding lengths `l0` and `l1`. The Green-Lagrange strain is then computed from the change in squared length:

```python
strain = 0.5 * (l1**2 - l0**2) / l0**2
```

That line is the key to the whole formulation. It is the large-deformation strain measure that makes the bar routines geometrically nonlinear.

`qebar` returns three pieces of information: the element force vector, the axial stress, and the strain. The force vector is built from the current element axis, not the initial axis, which is why the code uses `a1` in the final block:

```python
qe = (A * stress / l0) * np.vstack([-a1, a1])
```

The global wrappers `kbar` and `qbar` are vectorized. They collect all element nodes from the topology table, evaluate every element in one pass, and then scatter the element contributions back into the global arrays. The `np.einsum` calls are there to avoid explicit Python loops when forming repeated outer products like `a1 @ a1.T`. The `np.add.at` calls perform the scatter-add safely even when multiple elements write into the same global degrees of freedom.

The mass routines follow the same layout. `mebar` uses the compact material row `[A, E, rho]`, with `rho` defaulting to 1.0 when omitted. The `dof` argument controls the block size of the mass matrix, so the same kernel can represent 1D, 2D, or 3D bar-style nodes. If `lumped=True`, the total mass is split diagonally across the two nodes. Otherwise the consistent `[[2, 1], [1, 2]]` block is used.

Practical reading tip: when the comments mention "material stiffness" and "geometric stiffness", they mean the two additive pieces of the tangent matrix in `kebar`. When they mention "current element vector", they mean the element axis after deformation, not the original geometry.

## 2.3 Three-node triangles

The three-node triangle is the most compact 2D solid element in the library. It appears in both mechanical and scalar-field form:

- `ket3e`, `qet3e`, `kt3e`, `qt3e` for 2D elasticity.
- `ket3p`, `qet3p`, `kt3p`, `qt3p` for scalar potential or diffusion problems.
- `met3e`, `mt3e` for the triangle mass matrix.

The mechanical T3 element is the constant strain triangle. The implementation computes the edge-difference helpers once, then uses them to build a constant `B` matrix. Because the strain is constant over the element, no quadrature loop is needed. The geometry helper `_triangle_geometry` returns the edge-difference array and the area. The element stiffness then follows the compact formula `Ke = area * B.T @ D @ B`.

The comments in `ket3e` about `plane_strain` are important. The material row can optionally carry a type flag:

- `1` or omitted means plane stress.
- `2` means plane strain.

The code checks that flag directly and selects the constitutive matrix accordingly. That is why the docstring for `Ge` is short but the implementation is explicit about `props[2]`.

The internal-force kernel `qet3e` mirrors the stiffness kernel. It uses the same `B` matrix, multiplies it by the element displacement vector `Ue`, and then maps the strain to stress through `D`. The return value is a tuple of `(qe, stress, strain)`. The naming is consistent with the rest of the library, so `q` is the assembled force-like quantity, `S` is stress, and `E` is strain.

The global assembly functions `kt3e` and `qt3e` show the standard batched pattern used throughout the repository:

```python
nodes = topology[:, :3].astype(int) - 1
materials = as_float_array(G)[topology[:, -1].astype(int) - 1]
```

That pair of lines does two jobs. First, it removes the 1-based indexing inherited from the topology table. Second, it uses the material id stored in the last topology column to look up the correct row in `G`. After that, the code operates on batches of triangles using NumPy arrays of shape `(nel, 3, 2)`, `(nel, 3, 3)`, and `(nel, 6, 6)`.

The scalar triangle kernels `ket3p`, `qet3p`, `kt3p`, and `qt3p` reuse the same geometry but replace mechanics with diffusion or conduction. Here the material row is simpler: `Ge = [k]` or `Ge = [k, b]`. The optional second value is a reaction term. The comments about `B` and `D` still apply, but now `B` is a gradient operator and `D` is a conductivity matrix rather than an elastic constitutive matrix.

The mass functions `met3e` and `mt3e` deserve special attention because they show how the code stores thickness and density. The material row may carry extra columns:

- `Ge[3]` is thickness when present.
- `Ge[4]` is density when present.

If those columns are missing, the kernel defaults to 1.0. That means the mass code is permissive, but it also means you should be careful about the exact material-row format when building examples.

## 2.4 Quadrilateral elements

The quadrilateral family is where the code becomes more explicit about integration points, history arrays, and vectorization. The public routines are:

- `keq4e`, `qeq4e`, `kq4e`, `qq4e` for linear plane stress or plane strain elasticity.
- `keq4p`, `qeq4p`, `kq4p`, `qq4p` for scalar conductivity problems.
- `keq4eps`, `qeq4eps`, `kq4eps`, `qq4eps` for plane-stress elastoplasticity.
- `keq4epe`, `qeq4epe`, `kq4epe`, `qq4epe` for plane-strain elastoplasticity.
- `meq4e`, `mq4e` for mass matrices.

The shared geometry helpers `_q4_dN`, `_q4_B`, `_q4_gauss_points`, and `_q4_gp_index` tell you almost everything about the implementation style. `_q4_dN` evaluates parent-space derivatives, `_q4_gauss_points` returns the standard `2 x 2` Gauss rule, and `_q4_B` expands gradients into the mechanical strain-displacement layout. The Jacobian is never inverted explicitly. Instead, the code solves

```python
dN_global = np.linalg.solve(Jt, dN)
```

which is more stable and avoids building an intermediate inverse matrix.

`keq4e` is the cleanest example of the quadrilateral mechanics path. It loops over the four Gauss points, builds the element `B` matrix at each point, and accumulates `B.T @ D @ B * detJ * t`. The comments about thickness matter here because the element stiffness is scaled by `t` when thickness is present in the material row.

`qeq4e` uses the same integration loop but returns the local force vector, plus stress and strain values at each Gauss point. The `gp = _q4_gp_index(i, j)` helper is there to flatten the `2 x 2` tensor-product grid into a stable row order. This is why the history arrays in the plasticity kernels can be stored as one row per element and one block per Gauss point.

The scalar quadrilateral routines `keq4p`, `qeq4p`, `kq4p`, and `qq4p` follow the same structure but work on a single scalar unknown per node. The optional second material column is a reaction coefficient. When present, the code adds the `N.T @ (b * N)` contribution, which is the discrete form of the zero-order term.

The mass matrix `meq4e` is worth reading carefully because it shows both the consistent and lumped options in one place. The consistent path integrates `rho * t * N.T @ N * detJ` over the Gauss points. The lumped path first computes row sums, then rescales the diagonal so the total mass is preserved. That is a practical implementation choice rather than a theoretical one, and the comments in the source are pointing to exactly that tradeoff.

The plasticity kernels are the most stateful routines in the module.

- `keq4eps` and `qeq4eps` store plane-stress history data in arrays with width 16.
- `keq4epe` and `qeq4epe` store plane-strain history data in arrays with width 20.

The helper `_ensure_state_width` pads or reshapes those arrays so each element has enough room for every Gauss point. That is why the code passes `Se[i].reshape(4, 4)` or `Se[i].reshape(4, 5)` depending on the formulation. The extra columns are not random. They hold the evolving state variables needed by the return-mapping logic.

The `mtype` argument selects the yield model:

- `1` means von Mises.
- `2` means Drucker-Prager.

In the plane-stress path, `stressvm` and `stressdp` perform the stress correction after the trial elastic step. In the plane-strain path, the code uses a B-bar style modification to reduce volumetric locking. The important thing for a reader is not the exact algebraic formula in the docstring, but the flow: trial stress, yield test, corrected stress, updated history, then assembled force or tangent matrix.

Practical note on reading the code comments: when the source says "Gauss points" it means the four integration points used for the tensor-product rule. When it says "history tables", it means the per-element stored state that is passed back into the next nonlinear iteration.

## 2.5 Three-dimensional solids

The 3D solid family covers tetrahedra and hexahedra:

- `keT4e`, `qeT4e`, `kT4e`, `qT4e`, `meT4e`, `mT4e` for the 4-node tetrahedron.
- `keh8e`, `qeh8e`, `kh8e`, `qh8e`, `meh8e`, `mh8e` for the 8-node hexahedron.

The tetrahedral kernels are the most direct to read. The derivative matrix `dN_ref` is constant, so the Jacobian `J = dN_ref @ Xe` is evaluated once per element. That is why the T4 element is efficient: the shape-function gradients are constant, so there is no quadrature loop. The `B` matrix has shape `6 x 12`, reflecting the six Voigt strain components and the three displacement components at four nodes.

`keT4e` and `qeT4e` use the same `B` matrix and constitutive matrix `D`. `qeT4e` returns the element internal-force vector together with stress and strain vectors. `kT4e` and `qT4e` batch those operations across all tetrahedra by stacking the element coordinates into arrays of shape `(nel, 4, 3)`. The batch version uses `np.einsum` and broadcasted `np.linalg.solve` calls to keep the implementation fully vectorized.

The tetrahedral mass routines `meT4e` and `mT4e` are analytical. The consistent mass matrix uses the standard `rho * V / 20` block form, while the lumped version puts `rho * V / 4` on the diagonal. The code computes the volume from the determinant of the Jacobian, so the comments about `J` and `det(J)` directly match the implementation.

The hexahedral kernels are the 3D counterpart of the Q4 routines. They use `2 x 2 x 2` Gauss integration and a batched derivative helper `_hexa_dN_batch`. The flow is the same as in 2D, but the arrays are larger:

- `dN` has shape `(8, 3, 8)` for the eight Gauss points.
- `Jt` has shape `(8, 3, 3)` for one element, or `(nel, 8, 3, 3)` in the batched assembly path.
- `B` has shape `(8, 6, 24)` for one element, because there are six strain components and 24 displacement degrees of freedom.

`keh8e` and `qeh8e` are fully integrated element kernels. `kh8e` and `qh8e` perform global assembly. The code explicitly raises `NotImplementedError` for 20-node hexahedra, so the documented support is intentionally limited to the 8-node case.

The hexahedral mass routines `meh8e` and `mh8e` show the same design as the 2D mass kernels: consistent integration first, then an optional lumped path. The lumped path uses diagonal scaling to preserve total mass. That point is easy to miss if you only look at the function name, so it is worth checking the source comments when using the routine in a new example.

## 2.6 How to read the kernels

If you are trying to follow the source line by line, the recurring comment patterns are usually enough:

- "Compute shape function derivatives" means the code is building the parent-space gradient operator.
- "Jacobian" means the mapping from parent coordinates to physical coordinates.
- "Solve `J * dN_global = dN`" means the code is converting parent-space gradients into physical gradients without forming an explicit inverse.
- "Assemble into global matrix" means the function is using `assmk`, `assmq`, `element_dof_indices`, or `np.add.at` to scatter local contributions.
- "History array" means the kernel carries state from one nonlinear iteration to the next.

The naming conventions are consistent across the module:

- `k*` returns a stiffness-like matrix.
- `q*` returns an element force, flux, stress, or strain result.
- `m*` returns a mass matrix.
- `e` means element-level.
- `t3`, `q4`, `T4`, and `H8` identify the element family.

Once you learn that pattern, the element module becomes much easier to navigate. The code is not hiding complexity; it is laying out the mechanics in a repetitive form so that the same assembly logic can be reused across all element types.
