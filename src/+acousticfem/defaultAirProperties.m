function properties = defaultAirProperties()
%DEFAULTAIRPROPERTIES Air properties used by the original project examples.
properties = struct( ...
    'density', 1.23, ...                    % kg/m^3
    'dynamicViscosity', 1.95e-5, ...        % Pa s
    'thermalConductivity', 0.024, ...       % W/(m K)
    'heatCapacityRatio', 1.4, ...
    'ambientPressure', 1.013e5, ...         % Pa
    'specificHeatCapacity', 1004);           % J/(kg K)
end
