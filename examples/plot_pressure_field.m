function plot_pressure_field(result, plotTitle)
%PLOT_PRESSURE_FIELD Plot the real part of a nodal pressure field.
mesh = result.mesh;
figure('Name', plotTitle);
patch('Faces', mesh.elements, 'Vertices', mesh.nodes, ...
    'FaceVertexCData', real(result.pressure), 'FaceColor', 'interp', ...
    'EdgeColor', 'none');
axis equal tight;
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title(plotTitle, 'Interpreter', 'none');
end
