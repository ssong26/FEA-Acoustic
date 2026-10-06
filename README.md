# FEA-Acoustic

A small MATLAB finite-element solver for two-dimensional, frequency-domain acoustic pressure in homogeneous, porous, and heterogeneous media.

The author's original code is preserved in [archive/](archive/), under [archive/legacy/](archive/legacy/). In October 2026, the code was cleaned up and organized with assistance from Codex to make it easier to distribute. If you need the author's original implementation, please refer to the archived code. The original implementation dates to 2017; the maintained solver and examples now share one implementation.

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
