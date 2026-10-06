function result = graded_porous_material(makePlots, resolution)
%GRADED_POROUS_MATERIAL Air coupled to a continuously graded porous layer.
if nargin < 1, makePlots = true; end
if nargin < 2, resolution = [50, 50]; end
lengthX = 0.02;
heightY = 0.02;
frequency = 3550;
properties = acousticfem.defaultAirProperties();
mesh = acousticfem.rectangularMesh(lengthX, heightY, ...
    resolution(1), resolution(2));
centroids = (mesh.nodes(mesh.elements(:, 1), :) + ...
    mesh.nodes(mesh.elements(:, 2), :) + mesh.nodes(mesh.elements(:, 3), :)) / 3;
kAir = 2 * pi * frequency / 343;
elementK = repmat(kAir, size(mesh.elements, 1), 1);
inPorous = centroids(:, 1) > lengthX / 2;
radius = 1e-5 + 1e-3 * ...
    ((2 * centroids(:, 1) - lengthX) / lengthX) .* ...
    (centroids(:, 2) / heightY);
for e = find(inPorous).'
    elementK(e) = acousticfem.thermoviscousWavenumber( ...
        frequency, radius(e), properties);
end

tol = max(lengthX, heightY) * 1e-12;
fixed = find(abs(mesh.nodes(:, 1)) < tol);
result = acousticfem.solveHelmholtz(mesh, elementK, ...
    struct('nodeIndices', fixed, 'pressure', ones(size(fixed))));
result.elementWavenumbers = elementK;
if makePlots
    plot_pressure_field(result, 'Graded porous material');
end
end
