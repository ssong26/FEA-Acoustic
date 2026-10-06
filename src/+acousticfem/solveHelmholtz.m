function result = solveHelmholtz(mesh, elementWavenumbers, dirichlet)
%SOLVEHELMHOLTZ Solve the 2-D frequency-domain scalar Helmholtz equation.
%   dirichlet is a struct with nodeIndices and pressure fields. Unspecified
%   boundaries receive the natural homogeneous Neumann condition.

if ~isstruct(dirichlet) || ~isfield(dirichlet, 'nodeIndices') || ...
        ~isfield(dirichlet, 'pressure')
    error('acousticfem:InvalidBoundaryData', ...
        'dirichlet must contain nodeIndices and pressure.');
end
K = acousticfem.assembleSystem(mesh, elementWavenumbers);
nNodes = size(mesh.nodes, 1);
fixed = dirichlet.nodeIndices(:);
values = dirichlet.pressure(:);
if numel(fixed) ~= numel(values) || any(fixed < 1) || ...
        any(fixed > nNodes) || any(fixed ~= round(fixed)) || ...
        numel(unique(fixed)) ~= numel(fixed) || any(~isfinite(values))
    error('acousticfem:InvalidBoundaryData', ...
        'Dirichlet nodes and values must be finite, valid, and unique.');
end

free = setdiff((1:nNodes).', fixed);
p = complex(zeros(nNodes, 1));
p(fixed) = values;
if ~isempty(free)
    p(free) = K(free, free) \ (-K(free, fixed) * values);
end
result = struct('mesh', mesh, 'pressure', p, 'stiffness', K, ...
    'dirichletNodes', fixed, 'dirichletValues', values);
end
