function result = homogeneous_validation(makePlots, resolution)
%HOMOGENEOUS_VALIDATION Validate the solver against a 1-D analytic solution.
if nargin < 1, makePlots = true; end
if nargin < 2, resolution = [100, 10]; end

lengthX = 0.05;
heightY = 0.001;
frequency = 1000;
soundSpeed = 343;
k = 2 * pi * frequency / soundSpeed;
mesh = acousticfem.rectangularMesh(lengthX, heightY, resolution(1), resolution(2));
tol = max(lengthX, heightY) * 1e-12;
fixed = find(abs(mesh.nodes(:, 1) - lengthX) < tol);
solution = acousticfem.solveHelmholtz(mesh, k, ...
    struct('nodeIndices', fixed, 'pressure', ones(size(fixed))));

lineNodes = find(abs(mesh.nodes(:, 2) - heightY / 2) < tol);
[xLine, order] = sort(mesh.nodes(lineNodes, 1));
lineNodes = lineNodes(order);
numerical = solution.pressure(lineNodes);
analytical = cos(k * xLine) / cos(k * lengthX);
relativeError = norm(numerical - analytical) / max(norm(analytical), eps);
result = solution;
result.x = xLine;
result.numerical = numerical;
result.analytical = analytical;
result.relativeError = relativeError;
result.wavenumber = k;
fprintf('Homogeneous validation relative L2 error: %.6g\n', relativeError);

if makePlots
    figure('Name', 'Homogeneous validation');
    plot(xLine, real(numerical), '-', 'LineWidth', 1.5);
    hold on;
    plot(xLine, real(analytical), '--', 'LineWidth', 1.5);
    xlabel('x (m)'); ylabel('Pressure');
    legend('FEM', 'Analytical', 'Location', 'best');
    title('Homogeneous medium validation');
    grid on;
    plot_pressure_field(solution, 'Homogeneous pressure field');
end
end
