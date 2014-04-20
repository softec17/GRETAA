function Rk = R_Radar_altimeter(t,x,measurements,inputs,uncertainty)
% Calculates the measurement noise matrix for the radar altimeter
%
% Inputs
% t = current time
% x = current state vector
% measurements = a vector with the measurements at the current time
% inputs = a struct with information for state propagations
% uncertainty = a vector with uncertainty in the current measurements
%
% Output
% Rk = measurement sensitivity matrix (num_meas by num_meas)

Rk = diag(uncertainty).^2;

return