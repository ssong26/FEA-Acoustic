function result = thermoviscous_porous_layer(makePlots, resolution)
%THERMOVISCOUS_POROUS_LAYER Uniform thermoviscous porous equivalent-fluid layer.
if nargin < 1, makePlots = true; end
if nargin < 2, resolution = [100, 50]; end
lengthX = 0.1;
heightY = 0.05;
frequency = 1000;
poreRadius = 0.0005;
properties = acousticfem.defaultAirProperties();
k = acousticfem.thermoviscousWavenumber(frequency, poreRadius, properties);
mesh = acousticfem.rectangularMesh(lengthX, heightY, ...
    resolution(1), resolution(2));
tol = max(lengthX, heightY) * 1e-12;
fixed = find(abs(mesh.nodes(:, 1)) < tol);
result = acousticfem.solveHelmholtz(mesh, k, ...
    struct('nodeIndices', fixed, 'pressure', ones(size(fixed))));
result.wavenumber = k;
[~, effectiveDensity, effectiveBulkModulus] = ...
    acousticfem.thermoviscousWavenumber(frequency, poreRadius, properties);
result.effectiveDensity = effectiveDensity;
result.effectiveBulkModulus = effectiveBulkModulus;
if makePlots
    plot_pressure_field(result, 'Thermoviscous porous layer');
end
end
