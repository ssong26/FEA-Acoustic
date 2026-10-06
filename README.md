# FEA-Acoustic

A small MATLAB finite-element solver for two-dimensional, frequency-domain acoustic pressure in homogeneous, porous, and heterogeneous media.

The project grew out of an independent 2017 implementation. The maintained solver and examples now share one implementation; the original standalone files are preserved under `archive/legacy/` for provenance and are not part of the supported API.

## Scope

The current solver uses linear triangular (P1) elements for

$$
\nabla^2 p + k^2 p = 0,
$$

with the weak form

$$
\int_\Omega \nabla N^T\nabla p - k^2 N^T p\,d\Omega = 0.
$$

It supports prescribed pressure on selected nodes and the natural homogeneous Neumann condition on the remaining boundary. Porous examples use a thermoviscous effective-density/effective-bulk-modulus model. The circular-void example approximates the void by removing nearby mesh nodes and triangles; it is not an exact curved-boundary mesh.

## Requirements

- MATLAB R2021a or newer is recommended.
- The examples use MATLAB's built-in Delaunay triangulation, sparse matrices, and plotting functions. No additional toolbox is intentionally required.
- GNU Octave compatibility has not yet been verified.

## Run an example

From the repository root in MATLAB:

```matlab
addpath('src');
addpath('examples');
result = homogeneous_validation;
```

Other examples are `circular_void_scattering`, `thermoviscous_porous_layer`, `heterogeneous_materials`, and `graded_porous_material`. Each returns a struct containing the mesh, sparse stiffness matrix, and nodal complex pressure. Pass `false` as the first argument to suppress plots, and optionally pass `[nx, ny]` as the second argument to change mesh resolution.

## Tests

Run the lightweight solver and example checks from the repository root:

```matlab
addpath('tests');
run_tests;
```

The checks cover element-matrix properties, exact Dirichlet values, mesh-refinement behavior against the one-dimensional homogeneous solution, and finite solutions for each example. GitHub Actions runs the same checks in MATLAB on pushes and pull requests.

## Layout

```text
src/+acousticfem/   Mesh, element, material, assembly, and solve functions
examples/           Named, runnable physical scenarios
tests/               Small numerical and smoke checks
archive/legacy/      Original 2017 standalone scripts
archive/local/       Local-only archival files (ignored by Git)
```

## License and citation

The code is released under the MIT License. See [LICENSE](LICENSE). Citation metadata is in [CITATION.cff](CITATION.cff).
