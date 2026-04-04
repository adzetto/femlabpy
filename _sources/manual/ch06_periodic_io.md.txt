# Chapter 6: Periodic Boundaries & I/O

Advanced finite element analyses often require specialized boundary condition enforcement and integration with industrial meshing tools.

## 6.1 Periodic Boundary Conditions (Homogenization)

When simulating the mechanical response of a composite material, we typically analyze a Representative Volume Element (RVE). To ensure that the micro-scale deformations represent a continuous macro-scale material, Periodic Boundary Conditions (PBCs) must be applied. 

PBCs enforce that the deformation on opposite faces of the RVE are identical, offset only by the applied macroscopic strain tensor $\bar{\mathbf{\epsilon}}$. For two corresponding nodes $A^+$ (on the positive face) and $A^-$ (on the negative face), the displacement mapping is:

$$ \mathbf{u}^+ - \mathbf{u}^- = \bar{\mathbf{\epsilon}} \Delta \mathbf{x} $$

where $\Delta \mathbf{x} = \mathbf{x}^+ - \mathbf{x}^-$ is the physical distance vector between the paired nodes.

### Implementation in femlabpy

1. **Pairing Nodes:** The `find_periodic_pairs` function scans two arrays of node indices and links nodes that have matching coordinates along the boundary planes (within a numerical tolerance).
2. **Generating Constraints:** The `periodic_constraints` function converts these node pairs into linear constraint equations, assembling the $\mathbf{G}$ and $\mathbf{Q}$ matrices for the Lagrange multiplier solver.
3. **Homogenization Driver:** The `homogenize` function automates the entire process. It applies three independent macroscopic strain states (pure X-tension, pure Y-tension, and pure shear) to the RVE, solves the constrained equations using `solve_lag_general`, and computes the volume-averaged stress tensor $\langle \sigma \rangle$ for each state. The resultant vectors form the effective $3 \times 3$ elasticity matrix $\mathbf{C}_{eff}$.

## 6.2 Gmsh I/O Integration

`femlabpy` interfaces seamlessly with the open-source mesh generator [Gmsh](https://gmsh.info/). The `io.gmsh` module reads `.msh` files (versions 2.2 and 4.1) natively without requiring complex external dependencies.

The `load_gmsh2(filepath)` function parses the ASCII mesh block and returns a `GmshMesh` dataclass containing:
- `positions`: An $N \times 3$ array of spatial coordinates.
- `triangles`: An $E \times 4$ array of CST topologies.
- `quads`: An $E \times 5$ array of Q4 topologies.

This native array extraction bridges the gap between complex 3D CAD modeling and the flat-array `femlabpy` drivers instantly.