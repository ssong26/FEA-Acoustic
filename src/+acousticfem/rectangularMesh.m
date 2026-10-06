function mesh = rectangularMesh(lengthX, heightY, nx, ny, circularVoid)
%RECTANGULARMESH Delaunay mesh of a rectangle, optionally with a circular void.
%   nx and ny are the numbers of grid points in x and y. A void is a struct
%   with center=[x y] and radius fields; elements whose centroids lie inside
%   the circle are removed, matching the original example's approximation.

if nargin < 5
    circularVoid = [];
end
if any(~isfinite([lengthX, heightY])) || lengthX <= 0 || heightY <= 0 || ...
        nx < 2 || ny < 2 || nx ~= round(nx) || ny ~= round(ny)
    error('acousticfem:InvalidGrid', 'Lengths must be positive and nx, ny integers >= 2.');
end

[x, y] = meshgrid(linspace(0, lengthX, nx), linspace(0, heightY, ny));
nodes = [x(:), y(:)];
if ~isempty(circularVoid)
    if ~isfield(circularVoid, 'center') || ~isfield(circularVoid, 'radius') || ...
            numel(circularVoid.center) ~= 2 || circularVoid.radius <= 0
        error('acousticfem:InvalidVoid', 'Void needs a 2-D center and positive radius.');
    end
    d2 = sum((nodes - circularVoid.center(:).').^2, 2);
    nodes(d2 < circularVoid.radius^2, :) = [];
end
if size(nodes, 1) < 3
    error('acousticfem:InvalidGrid', 'Mesh has fewer than three nodes.');
end

elements = delaunay(nodes(:, 1), nodes(:, 2));
if ~isempty(circularVoid)
    centroids = (nodes(elements(:, 1), :) + nodes(elements(:, 2), :) + ...
        nodes(elements(:, 3), :)) / 3;
    d2 = sum((centroids - circularVoid.center(:).').^2, 2);
    elements(d2 < circularVoid.radius^2, :) = [];
end
if isempty(elements)
    error('acousticfem:InvalidGrid', 'Mesh contains no elements.');
end
mesh = struct('nodes', nodes, 'elements', elements);
end
