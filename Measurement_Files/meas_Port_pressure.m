function [GaugePress] = meas_Port_pressure(t,x,PressureInputs,inputs)
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
lat = x(2);
lon = x(3);
u = x(4); % Velocity in North direction (of NED frame)
v = x(5); % Velocity in East direction (of NED frame)
w = x(6); % Velocity in Down direction (of NED frame)
% Quaternions
q_J2000_2_DS_0 = x(7);
q_J2000_2_DS_1 = x(8);
q_J2000_2_DS_2 = x(9);
q_J2000_2_DS_3 = x(10);
q_J2000_2_MCMF_0 = x(11);
q_J2000_2_MCMF_1 = x(12);
q_J2000_2_MCMF_2 = x(13);
q_J2000_2_MCMF_3 = x(14);
% Atmosphere
rhoinf = x(15);
Pinf = x(16);

alt = r - inputs.Re;
Vinf = norm([u;v;w]);

R_J2000_2_DS = quat2dcm([q_J2000_2_DS_0,q_J2000_2_DS_1,q_J2000_2_DS_2,q_J2000_2_DS_3]); R_DS_2_J2000 = R_J2000_2_DS';
R_J2000_2_MCMF = quat2dcm([q_J2000_2_MCMF_0,q_J2000_2_MCMF_1,q_J2000_2_MCMF_2,q_J2000_2_MCMF_3]);
R_MCMF_2_LH = R2(3*pi/2-lat)*R3(lon);
% Rotation vector from veh. to body
R_DS_LH = R_MCMF_2_LH*R_J2000_2_MCMF*R_DS_2_J2000; R_LH_DS = R_DS_LH';
R_DS_2_Body = [0 0 1;0 1 0;-1 0 0];
Rveh2body = R_DS_2_Body*R_LH_DS;

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