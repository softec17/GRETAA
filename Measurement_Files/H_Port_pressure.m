function [Hk,zk] = H_Port_pressure(t,x,measurements,inputs)
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


PressureInputs = inputs.PressureInputs;
[PortPressnom] = meas_Port_pressure(t,x,PressureInputs,inputs);
PortPressnom = PortPressnom(:);
num_x = length(x);
num_meas = length(measurements);

% Preallocate space
Hk = zeros(num_meas,num_x);

% Numerical differentiation step size
delh = 1e-5;

parfor jj = 15:16
    xplus = x;
    xminus = x;
    delx = delh*abs(x(jj)) + eps;
    xplus(jj) = xplus(jj) + delx;
    xminus(jj) = xminus(jj) - delx;
    [PortPressplus] = meas_Port_pressure(t,xplus,PressureInputs,inputs);
    [PortPressminus] = meas_Port_pressure(t,xminus,PressureInputs,inputs);
    Hk(:,jj) = (PortPressplus'-PortPressminus')/(2*delx);
end

measurements = measurements(:);
zk = measurements - PortPressnom;

return