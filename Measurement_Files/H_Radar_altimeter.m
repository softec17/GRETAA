function [Hk,zk] = H_Radar_altimeter(t,x,measurements,inputs)
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

% Radius (meters)
Radius = x(1);
u = x(4);
v = x(5);
w = x(6);
% Altitude (meters)
h = Radius - inputs.Re;
V = norm([u;v;w]);

% Altimeter
zk = measurements - h;
% Length of states
num_x = length(x);
% Preallocate space
Hk = zeros(num_meas,num_x);
% Allocate measurement sensitivity matrix
Hk(1,1) = 1;

% % Slant range
% gamma = atan(-w/sqrt(u^2+v^2));
% slantrange = h/sin(gamma);
% zk = measurements - slantrange;
% % Length of states
% num_x = length(x);
% % Preallocate space
% Hk = zeros(num_meas,num_x);
% % Allocate measurement sensitivity matrix
% Hk(1,1) = sin(gamma);
% Hk(1,4) = (h/w)*(V^-0.5)*(u);
% Hk(1,5) = (h/w)*(V^-0.5)*(v);
% Hk(1,6) = h/V - (h*V)/w^2;

return