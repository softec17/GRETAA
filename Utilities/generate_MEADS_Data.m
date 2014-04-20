function [Pressure_measurement] = generate_MEADS_Data(actual,PressureInputs)
% This will be static pressure measurement (not gauge pressure)
% Inputs
% actual = simulated POST2 output struct
% PressureInputs = struct with pressure port information
% Output
% Pressure_measurement = vector of measurements for time history

% Number of ports
N_port = length(PressureInputs.clock);

% Preallocating space
Pressure_measurement = zeros(length(actual.time),N_port);

% Data table containing CFD data
aerodata = PressureInputs.aero_data;
% Clock and cone angles (rad)
clock = PressureInputs.clock;
cone = PressureInputs.cone;
% Actual state variables from the POST2 data
P_inf = actual.pres;
rho_inf = actual.dens;
mach_inf = actual.mach;
alpha_inf = actual.alpha; % (deg)
beta_inf = actual.beta;   % (deg)
% soundspeed = actual.cs;   % Speed of sound (m/s)  
velr = actual.velr;

disp('Generating MEADS dataset');
parfor ii = 1:length(actual.time)
    % Angle of attack and sideslip angle
    alpha = deg2rad(alpha_inf(ii));
    beta = deg2rad(beta_inf(ii));
    % Freestream Mach number
    mach = mach_inf(ii);
%     speedSound = soundspeed(ii);
    % Freestream velocity
%     Vinf = mach*speedSound;
    Vinf = velr(ii);
    
    % Look up Cp from CFD table
    % Cp is an [num_ports,1] array
    [Cp] = CFDInterpolate(clock,cone,alpha,beta,mach,aerodata);
    GaugePress = Cp.*(0.5).*rho_inf(ii).*Vinf.^2;
%     PortPress = GaugePress + P_inf(ii);
    PortPress = GaugePress;
    Pressure_measurement(ii,:) = PortPress;
end

return