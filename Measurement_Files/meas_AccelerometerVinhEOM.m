function [accel] = meas_AccelerometerVinhEOM(t,x,inputs)
% Inputs
% t = time (sec)
% x = states
% inputs = structure containing information for state propagation
%
% Output
% range = range from station to s/c (m)
% azimuth = heading angle from station to s/c (deg)
% elevation = elevation angle from station to s/c (deg)

%% Current states
V = x(4);
gamma      = x(5);  %rad, relative flight path angle
psi        = x(6);  %rad, relative heading angle
q0 = x(7);
q1 = x(8);
q2 = x(9);
q3 = x(10);
% rhoinf = x(12);
rhoinf = interp1q(inputs.rhoTime,inputs.rhoHist,t);

% Other secondary states
u = V*cos(gamma)*sin(psi);
v = V*cos(gamma)*cos(psi);
w = -V*sin(gamma);

%% Vehicle properties
% Local speed of sound (m/s)
speedSound = interp1q(inputs.actual.time,inputs.actual.cs,t);

% Calculate Mach infinity (or Mach vehicle when W = 0)
mach = V/speedSound;

% Rotation vector from veh. to body
% Rveh2body = R_veh2body_wq(q0,q1,q2,q3);
Rveh2body = quat2dcm([q0,q1,q2,q3]);

% Vehicle velocity wrt ground (and wind when W = 0) in body frame
vel_body = Rveh2body*[u;v;w];

% Calculate angle-of-attack (deg) and sideslip angle (deg)
[alpha,beta] = Vw_2_alphabeta(vel_body);
alpha = deg2rad(alpha);
beta = deg2rad(beta);

% [CA,~,~,~,~] = aero_interp2(inputs.AOA_grid,inputs.Mach_grid,inputs.data_grid,alpha,beta,mach);
% CA = interp1q(inputs.CDtime,inputs.CA,t);

A = inputs.RefArea;

if length(inputs.mass) ~= 1
    mass = interp1q(inputs.CDtime,inputs.mass,t);
else
    mass = inputs.mass;
end

% accel(1) = 0.5*rhoinf*V^2*CA*A/mass;

% rhoinf = interp1q(inputs.actual.time,inputs.actual.dens,t);
% V = interp1q(inputs.actual.time,inputs.actual.velr,t);
% alpha = interp1q(inputs.actual.time,inputs.actual.alpha,t)*pi/180;
% beta = interp1q(inputs.actual.time,inputs.actual.beta,t)*pi/180;
% mach = interp1q(inputs.actual.time,inputs.actual.mach,t);
[CA] = aero_interp2_CA(inputs.AOA_grid,inputs.Mach_grid,inputs.data_grid,alpha,beta,mach);
accel(1) = 0.5*rhoinf*V^2*CA*A/mass;


return