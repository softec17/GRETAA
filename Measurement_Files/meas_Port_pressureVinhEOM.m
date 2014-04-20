function [GaugePress] = meas_Port_pressureVinhEOM(t,x,PressureInputs,inputs)
% Creates the pressure measurements at the forebody pressure transducers
% Inputs
% t = current time (not used in the calculations)
% x = current estimated state (see EOM file for the full list)
% PressureInputs = struct containing the look-up table with Cp data
% Outputs
% PortPress = stagnation pressure at port
% GaugePress = dynamic pressure at port

% Current states
r = x(1);
V = x(4);
gamma      = x(5);  %rad, relative flight path angle
psi        = x(6);  %rad, relative heading angle
% Quaternions
q0 = x(7);
q1 = x(8);
q2 = x(9);
q3 = x(10);
% Atmosphere
rhoinf = x(12);

% Other secondary states
u = V*cos(gamma)*sin(psi);
v = V*cos(gamma)*cos(psi);
w = -V*sin(gamma);
alt = r - inputs.Re;
Vinf = V;

% Rotation matrix from vehicle-carried horizontal frame to body frame
Rveh2body = quat2dcm([q0,q1,q2,q3]);

% Vehicle velocity wrt ground (and wind when W = 0) in body frame
vel_body = Rveh2body*[u;v;w];
% Calculate angle-of-attack (deg) and sideslip angle (deg)
[alpha,beta] = Vw_2_alphabeta(vel_body);
alpha = deg2rad(alpha);
beta = deg2rad(beta);
% Clock and cone angles (rad)
clock = PressureInputs.clock;
cone = PressureInputs.cone;

% Local speed of sound (m/s)
SoundSpeed = interp1(inputs.PressureInputs.speedSound.altlist,inputs.PressureInputs.speedSound.speedsoundlist,alt);

% Calculate Mach infinity (or Mach vehicle when W = 0)
mach = Vinf/SoundSpeed;
% Data table containing CFD data
aerodata = PressureInputs.aero_data;

% Look up Cp from CFD table
% Cp is an [num_ports,1] array
% alpha = interp1q(inputs.actual.time,inputs.actual.alpha,t).*pi/180;
% beta = interp1q(inputs.actual.time,inputs.actual.beta,t).*pi/180;
% mach = interp1q(inputs.actual.time,inputs.actual.mach,t);
[Cp] = CFDInterpolate(clock,cone,alpha,beta,mach,aerodata);
GaugePress = Cp.*(0.5).*rhoinf.*Vinf.^2;

return