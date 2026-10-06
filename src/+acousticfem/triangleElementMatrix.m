function Ke = triangleElementMatrix(xy, wavenumber)
%TRIANGLEELEMENTMATRIX Exact P1 Helmholtz element matrix on a linear triangle.
%   The weak form is integral(grad(N)'*grad(N) - k^2*N'*N) dOmega.

if ~isequal(size(xy), [3, 2]) || any(~isfinite(xy(:)))
    error('acousticfem:InvalidTriangle', 'xy must be a finite 3-by-2 array.');
end
if ~isscalar(wavenumber) || ~isfinite(wavenumber)
    error('acousticfem:InvalidWavenumber', 'wavenumber must be a finite scalar.');
end

x = xy(:, 1);
y = xy(:, 2);
detJ = (x(2) - x(1)) * (y(3) - y(1)) - ...
       (x(3) - x(1)) * (y(2) - y(1));
area = abs(detJ) / 2;
if area <= eps(max(1, max(abs(xy(:)))))^2
    error('acousticfem:DegenerateTriangle', 'Triangle area is zero or numerically negligible.');
end

% Each row is the constant physical gradient of one shape function.
gradN = [y(2) - y(3), x(3) - x(2); ...
         y(3) - y(1), x(1) - x(3); ...
         y(1) - y(2), x(2) - x(1)] / (2 * area);
mass = (area / 12) * (ones(3) + eye(3));
Ke = area * (gradN * gradN.') - (wavenumber^2) * mass;
end
