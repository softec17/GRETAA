function xdot = SixDOFEntryEOMVinhEOMMSLEOMAccel(t,x,inputs)
% Equations of motions for 6-DOF where position is r, lat, lon and
% velocity is v, gamma(FPA) and psi (heading angle) from the local
% horizontal frame

% disp(t)

%% Define state variables
r          = x(1);  %m,   planetocentric radius
lat   = x(2);  %rad, geocentric lat
lon = x(3);    %rad, longitude
u = x(4);   % North velocity inertial
v = x(5);   % East velocity inertial
w = x(6);   % Down velocity inertial
qJDS0 = x(7);          % Quaternions J2000_DS
qJDS1 = x(8);
qJDS2 = x(9);
qJDS3 = x(10);
qJMCMF0 = x(11);
qJMCMF1 = x(12);
qJMCMF2 = x(13);
qJMCMF3 = x(14);
rho = x(15);
Pinf = x(16);
alt = r - inputs.Re;

q_J2000_DS = [qJDS0,qJDS1,qJDS2,qJDS3];
q_J2000_MCMF = [qJMCMF0,qJMCMF1,qJMCMF2,qJMCMF3];
R_J_DS = quat2dcm(q_J2000_DS); R_DS_2_J = R_J_DS';
R_J_MCMF = quat2dcm(q_J2000_MCMF);
R_MCMF_LH = R2(3*pi/2-lat)*R3(lon);
R_DS_LH = R_MCMF_LH*R_J_MCMF*R_DS_2_J; R_LH_DS = R_DS_LH';
q_LH_2_DS = dcm2quat(R_LH_DS);

%% Accelerations
% Acceleration
acc_list = inputs.accel;
acc_time_list = inputs.time_accel;
acc = interp1q(acc_time_list,acc_list,t); acc = acc(:);

% IMU converted to LH frame
acc = R_DS_LH*acc;

% Gyro
gyro = interp1q(inputs.time_gyro,inputs.gyro,t); gyro = gyro(:);    % In DS frame

% Planet rotation rate
omega_planet = inputs.omegaplanet;

% Gravity terms
mu = inputs.mu;
J2 = inputs.J2;

%% Equations of motion
% Derivative terms
r_dot = -w;
lat_dot = u/r;
lon_dot = v/(r*cos(lat)) - omega_planet;
udot = acc(1) - (3*mu*J2/(2*r^4))*sin(2*lat) + (u*w-v^2*tan(lat))/r ;
vdot = acc(2) + (u*v*tan(lat)+v*w)/r ;     %- omega_planet*r*cos(lat) % Planet relative velocity
wdot = acc(3) + mu/r^2 - (3*mu*J2/(2*r^4))*(2-3*cos(lat)^2) - (u^2 + v^2)/r ;
q_J2000_2_DS_dot = quatEOM_general(t,q_J2000_DS,gyro);
q_J2000_2_MCMF_dot = quatEOM(t,q_J2000_MCMF,omega_planet);
g_local_r = mu/r^2;
rhoinfdot = (rho^2)*g_local_r*w/Pinf;
Pinfdot = rho*g_local_r*w;

xdot = [r_dot;lat_dot;lon_dot;udot;vdot;wdot;q_J2000_2_DS_dot;q_J2000_2_MCMF_dot;rhoinfdot;Pinfdot];

if sum(isnan(xdot))>=1
    keyboard
end

return
