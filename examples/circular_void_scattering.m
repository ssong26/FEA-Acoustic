function result = circular_void_scattering(makePlots, resolution)
%CIRCULAR_VOID_SCATTERING Acoustic field around a circular void in air.
if nargin < 1, makePlots = true; end
if nargin < 2, resolution = [100, 50]; end
lengthX = 0.1;
heightY = 0.05;
frequency = 10000;
soundSpeed = 343;
void = struct('center', [lengthX / 2, heightY / 2], ...
    'radius', heightY / 4);
mesh = acousticfem.rectangularMesh(lengthX, heightY, ...
    resolution(1), resolution(2), void);
tol = max(lengthX, heightY) * 1e-12;
fixed = find(abs(mesh.nodes(:, 1)) < tol);
result = acousticfem.solveHelmholtz(mesh, 2 * pi * frequency / soundSpeed, ...
    struct('nodeIndices', fixed, 'pressure', ones(size(fixed))));
result.void = void;
if makePlots
    plot_pressure_field(result, 'Circular void scattering');
end
end
