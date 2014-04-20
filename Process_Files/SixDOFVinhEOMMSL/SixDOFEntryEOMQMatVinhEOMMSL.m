function QMat = SixDOFEntryEOMQMatVinhEOMMSL(t,x,inputs)

%% Initialize matrices
n_state = length(x);
QMat = zeros(n_state);

%% Current states

% Quaternions
q0 = x(7);
q1 = x(8);
q2 = x(9);
q3 = x(10);

% Atmospheric states
PressInf = x(11);
DensInf = x(12);

%% Meas. uncer -> process noise IMU angular rate uncertainties
gyroU = inputs.IMU.gyroUncertainty;

%% State uncer -> process noise
% Velocity component process noise
% R_veh2body = R_veh2body_wq(q0,q1,q2,q3);
% R_body2veh = R_veh2body';
% QMat(4:6,4:6) = diag(R_body2veh*Qvec(1:3));
velNoise = inputs.Process.vel;
gamNoise = inputs.Process.gam;
psiNoise = inputs.Process.psi;
QMat(4:6,4:6) = diag([velNoise;gamNoise;psiNoise]);

%dqdot/dgyro
QMat(7:10,7:10) = diag(0.5.*qmult_left(q0,q1,q2,q3)*[0;gyroU]);

% dPressdot/dnoise and dDensdot/dnoise
% Pressure and density process noise is a fraction of their current values
PresUMult = inputs.Process.press;
DensUMult = inputs.Process.dens;
if ((PressInf/1000)>0.1)
    PresUMult = PresUMult*1e1;
    DensUMult = DensUMult*1e1;
end
PresU = PresUMult*PressInf;
DensU = DensUMult*DensInf;
% PresU = PresUMult;
% DensU = DensUMult;
QMat(11:12,11:12) = diag([PresU;DensU]);

% Position noise
radiusNoise = inputs.Process.rad;
latNoise = inputs.Process.lat;
lonNoise = inputs.Process.lon;
QMat(1:3,1:3) = diag([radiusNoise;latNoise;lonNoise]);

%% Turn 1 sigma uncertainty to process noise covariance
QMat = QMat.^2;

return