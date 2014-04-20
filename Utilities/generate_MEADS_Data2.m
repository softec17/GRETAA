function [Pressure_measurement] = generate_MEADS_Data2(actual,PressureInputs,indices)
% This will be static pressure measurement (not gauge pressure)
% Inputs
% actual = simulated POST2 output struct
% PressureInputs = struct with pressure port information
% Output
% Pressure_measurement = vector of measurements for time history

% Number of ports
N_port = length(PressureInputs.clock);

% Preallocating space
Pressure_measurement = zeros(length(actual.time(indices)),N_port);

% Data table containing CFD data
aerodata = PressureInputs.aero_data;
% Clock and cone angles (rad)
clock = PressureInputs.clock;
cone = PressureInputs.cone;
% Actual state variables from the POST2 data
P_inf = actual.pres(indices);
rho_inf = actual.dens(indices);
mach_inf = actual.mach(indices);
alpha_inf = actual.alpha(indices); % (deg)
beta_inf = actual.beta(indices);   % (deg)
% soundspeed = actual.cs;   % Speed of sound (m/s)
velr = actual.velr(indices);

% Angle of attack and sideslip angle
alpha = deg2rad(alpha_inf);
beta = deg2rad(beta_inf);
% Freestream Mach number
mach = mach_inf;
Vinf = velr;

parfor ii = 1:length(indices)
    % Look up Cp from CFD table
    % Cp is an [num_ports,1] array
    [Cp] = CFDInterpolate(clock,cone,alpha(ii),beta(ii),mach(ii),aerodata);
    GaugePress = Cp.*(0.5).*rho_inf(ii).*Vinf(ii).^2;
    PortPress = GaugePress + P_inf(ii);
    Pressure_measurement(ii,:) = PortPress;
end

return