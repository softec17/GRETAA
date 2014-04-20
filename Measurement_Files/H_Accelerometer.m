function [Hk,zk] = H_Accelerometer(t,x,measurements,inputs)
% Calculates the measurement sensitivity matrix and the residual in
% predicted and actual measurements
%
% Inputs
% t = current time
% x = current state vector (this function doesn't care about the size)
% measurements = vector with the measurements
% inputs = struct containing information needed for state propagation and
%          EKF update. Contents depend on what application is being run
%
% Outputs
% Hk = measurement sensisitivity matrix (n_meas by n_states)
% zk = measurement residual (btw. predicted and actual measurements)

% Number of measurements
num_meas = length(measurements);

% Length of states
num_x = length(x);

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
Rveh2body = R_veh2body_wq(q0,q1,q2,q3);

% Vehicle velocity wrt ground (and wind when W = 0) in body frame
vel_body = Rveh2body*[u;v;w];

% Calculate angle-of-attack (deg) and sideslip angle (deg)
[alpha,beta] = Vw_2_alphabeta(vel_body);
alpha = deg2rad(alpha);
beta = deg2rad(beta);

[CA] = aero_interp2_CA(inputs.AOA_grid,inputs.Mach_grid,inputs.data_grid,alpha,beta,mach);
% CA = interp1q(inputs.CDtime,inputs.CA,t);

A = inputs.RefArea;

if length(inputs.mass) ~= 1
    mass = interp1q(inputs.CDtime,inputs.mass,t);
else
    mass = inputs.mass;
end

%% Accelerations (accel_x,accel_y,accel_z) => (-CA,-CY,-CN)
% accel(1) = 0.5*rhoinf*V^2*CA*A/mass;

% rhoinf = interp1q(inputs.actual.time,inputs.actual.dens,t);
% V = interp1q(inputs.actual.time,inputs.actual.velr,t);
% alpha = interp1q(inputs.actual.time,inputs.actual.alpha,t)*pi/180;
% beta = interp1q(inputs.actual.time,inputs.actual.beta,t)*pi/180;
% mach = interp1q(inputs.actual.time,inputs.actual.mach,t);
% [CA,~,~,~,~] = aero_interp2(inputs.AOA_grid,inputs.Mach_grid,inputs.data_grid,alpha,beta,mach);
[CA] = aero_interp2_CA(inputs.AOA_grid,inputs.Mach_grid,inputs.data_grid,alpha,beta,mach);
accel(1) = 0.5*rhoinf*V^2*CA*A/mass;

% Measurement residual
zk = measurements(:) - accel(:);

%% Jacobian
% Preallocate space
Hk = zeros(1,num_x);
% Hk(:,4) = rhoinf*V*CA*A/mass;

% delh = 1e-11;
% parfor jj = 1:num_x
%     xplus = x; xminus = x;
%     dx = abs(x(jj))*delh + eps;
%     xplus(jj) = xplus(jj) + dx;
%     xminus(jj) = xminus(jj) - dx;
%     xdotplus = meas_Accelerometer(t,xplus,inputs);
%     xdotminus = meas_Accelerometer(t,xminus,inputs);
%     Hk(:,jj) = (xdotplus-xdotminus)/(2*dx);
% end
% Hk(:,1) = (0.5*V*CA*A/mass)*rhoinf*(1/11029.01);
Hk(:,4) = rhoinf*V*CA*A/mass;
Hk(:,12) = 0.5*V^2*CA*A/mass;

% Hk = 1e-1.*Hk;

return