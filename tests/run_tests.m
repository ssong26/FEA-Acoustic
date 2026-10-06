function run_tests()
%RUN_TESTS Lightweight MATLAB checks for assembly, boundary data, and examples.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'examples'));
test_element_matrix();
test_boundary_conditions_and_convergence();
test_remaining_examples();
fprintf('All FEA-Acoustic checks passed.\n');
end

function test_element_matrix()
xy = [0, 0; 1, 0; 0, 1];
K0 = acousticfem.triangleElementMatrix(xy, 0);
assert(norm(K0 - K0.', 'fro') < 1e-12, 'Element matrix must be symmetric.');
assert(norm(K0 * ones(3, 1)) < 1e-12, 'Zero-wavenumber element must annihilate constants.');
Kc = acousticfem.triangleElementMatrix(xy, 2 + 0.5i);
assert(all(isfinite(Kc(:))), 'Complex element matrix must be finite.');
end

function test_boundary_conditions_and_convergence()
coarse = homogeneous_validation(false, [32, 6]);
fine = homogeneous_validation(false, [64, 12]);
assert(all(isfinite(fine.pressure)), 'Homogeneous solution must be finite.');
assert(fine.relativeError < coarse.relativeError, ...
    'Refining the mesh should reduce the homogeneous validation error.');
assert(max(abs(fine.pressure(fine.dirichletNodes) - ...
    fine.dirichletValues)) < 1e-12, 'Prescribed pressure must be imposed exactly.');
end

function test_remaining_examples()
results = {circular_void_scattering(false, [24, 12]), ...
    thermoviscous_porous_layer(false, [20, 10]), ...
    heterogeneous_materials(false, [20, 14]), ...
    graded_porous_material(false, [20, 14])};
for k = 1:numel(results)
    assert(all(isfinite(results{k}.pressure)), ...
        'Each example must return a finite pressure field.');
end
end
