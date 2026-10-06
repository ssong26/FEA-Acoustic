function [k, effectiveDensity, effectiveBulkModulus] = ...
        thermoviscousWavenumber(frequency, poreRadius, properties)
%THERMOVISCOUSWAVENUMBER Effective-fluid coefficient from the original slit model.
if nargin < 3
    properties = acousticfem.defaultAirProperties();
end
if ~isscalar(frequency) || ~isfinite(frequency) || frequency <= 0 || ...
        ~isscalar(poreRadius) || ~isfinite(poreRadius) || poreRadius <= 0
    error('acousticfem:InvalidMaterialData', ...
        'Frequency and pore radius must be positive finite scalars.');
end

rho = properties.density;
mu = properties.dynamicViscosity;
cp = properties.specificHeatCapacity;
conductivity = properties.thermalConductivity;
gamma = properties.heatCapacityRatio;
p0 = properties.ambientPressure;
Pr = mu * cp / conductivity;
omega = 2 * pi * frequency;
s = sqrt(omega * rho * poreRadius^2 / mu);
rootI = sqrt(1i);
effectiveDensity = rho / (1 - tanh(s * rootI) / (s * rootI));
effectiveBulkModulus = gamma * p0 / (1 + (gamma - 1) * ...
    tanh(sqrt(Pr) * s * rootI) / (sqrt(Pr) * s * rootI));
k = omega * sqrt(effectiveDensity / effectiveBulkModulus);
end
