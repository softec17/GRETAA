function Rk = R_Port_pressure(t,x,measurements,inputs,uncertainty)
% Calculates the measurement noise matrix for the port pressure measurement
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

% Current states
r = x(1);
u = x(4); % Velocity in North direction (of NED frame)
v = x(5); % Velocity in North direction (of NED frame)
w = x(6); % Velocity in North direction (of NED frame)

alt = r - inputs.Re;
V = norm([u;v;w]);

% Freestream pressure
% PInf = x(11);
% % Freestream density
% DensInf = x(12);
% DensInf = interp1q(inputs.rhoTime,inputs.rhoHist,t);

% Other parameters
SoundSpeed = interp1(inputs.PressureInputs.speedSound.altlist,inputs.PressureInputs.speedSound.speedsoundlist,alt);
Mach = V/SoundSpeed;

% Pressure measurements
PortPress = measurements;

% Pressure measurement uncertainty
% Uncertainty is a multiplier to the current pressure value
% Will use mean port pressure as the "current" pressure value
MeanPortPress = mean(PortPress);
PortPressUncertainty = uncertainty.*MeanPortPress;

%     Calculates the measurement uncertainty covariance
if Mach >= 1.5
    Rk = diag(PortPressUncertainty).^2;
else  % Ignore data after Mach 1.5
    uncertainty(:) = uncertainty.*1000;
    PortPressUncertainty = uncertainty.*MeanPortPress;
    Rk = diag(PortPressUncertainty).^2;
end

return