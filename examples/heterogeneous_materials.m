function result = heterogeneous_materials(makePlots, resolution)
%HETEROGENEOUS_MATERIALS Air and two porous materials in a piecewise domain.
if nargin < 1, makePlots = true; end
if nargin < 2, resolution = [50, 50]; end
lengthX = 0.02;
heightY = 0.02;
frequency = 3550;
properties = acousticfem.defaultAirProperties();
mesh = acousticfem.rectangularMesh(lengthX, heightY, ...
    resolution(1), resolution(2));

kAir = 2 * pi * frequency / 343;
kPorous1 = acousticfem.thermoviscousWavenumber(frequency, 0.00002, properties);
kPorous2 = acousticfem.thermoviscousWavenumber(frequency, 0.01, properties);
centroids = (mesh.nodes(mesh.elements(:, 1), :) + ...
    mesh.nodes(mesh.elements(:, 2), :) + mesh.nodes(mesh.elements(:, 3), :)) / 3;
elementK = repmat(kAir, size(mesh.elements, 1), 1);
inPorous = centroids(:, 1) > lengthX / 2;
inMiddleBand = centroids(:, 2) > heightY / 3 & ...
    centroids(:, 2) < 2 * heightY / 3;
elementK(inPorous & ~inMiddleBand) = kPorous1;
elementK(inPorous & inMiddleBand) = kPorous2;

tol = max(lengthX, heightY) * 1e-12;
fixed = find(abs(mesh.nodes(:, 1)) < tol);
result = acousticfem.solveHelmholtz(mesh, elementK, ...
    struct('nodeIndices', fixed, 'pressure', ones(size(fixed))));
result.elementWavenumbers = elementK;
if makePlots
    plot_pressure_field(result, 'Piecewise heterogeneous materials');
end
end
