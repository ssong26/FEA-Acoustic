function K = assembleSystem(mesh, elementWavenumbers)
%ASSEMBLESYSTEM Assemble a sparse P1 Helmholtz matrix from triangular elements.

if ~isstruct(mesh) || ~isfield(mesh, 'nodes') || ~isfield(mesh, 'elements')
    error('acousticfem:InvalidMesh', 'mesh must contain nodes and elements.');
end
nodes = mesh.nodes;
elements = mesh.elements;
if size(nodes, 2) ~= 2 || size(elements, 2) ~= 3 || ...
        any(elements(:) < 1) || any(elements(:) > size(nodes, 1))
    error('acousticfem:InvalidMesh', 'Expected 2-D nodes and valid 3-node elements.');
end

nElements = size(elements, 1);
if isscalar(elementWavenumbers)
    elementWavenumbers = repmat(elementWavenumbers, nElements, 1);
else
    elementWavenumbers = elementWavenumbers(:);
end
if numel(elementWavenumbers) ~= nElements || any(~isfinite(elementWavenumbers))
    error('acousticfem:InvalidMaterialData', ...
        'Provide one finite wavenumber or one per element.');
end

nNodes = size(nodes, 1);
I = zeros(9 * nElements, 1);
J = zeros(9 * nElements, 1);
V = zeros(9 * nElements, 1);
for e = 1:nElements
    conn = elements(e, :);
    Ke = acousticfem.triangleElementMatrix(nodes(conn, :), elementWavenumbers(e));
    [ii, jj] = ndgrid(conn, conn);
    offset = 9 * (e - 1) + (1:9);
    I(offset) = ii(:);
    J(offset) = jj(:);
    V(offset) = Ke(:);
end
K = sparse(I, J, V, nNodes, nNodes);
end
