function xdot = SixDOFEntryEOMVinhEOMAdaptive(t,x,inputs)
% Equations of motions for 6-DOF where position is r, lat, lon and
% velocity is v, gamma(FPA) and psi (heading angle) from the local
% horizontal frame

%% Define state variables
r          = x(1);  %m,   planetocentric radius
lat   = x(2);  %rad, geocentric lat
V_relative = x(4);  %m/s, relative velocity
gamma      = x(5);  %rad, relative flight path angle
psi        = x(6);  %rad, relative heading angle
q0 = x(7);          % Quaternions
q1 = x(8);
q2 = x(9);
q3 = x(10);
Pinf = x(11);
% rho = x(12);
rho = interp1q(inputs.rhoTime,inputs.rhoHist,t);

% Other secondary states
u = V_relative*cos(gamma)*sin(psi);
v = V_relative*cos(gamma)*cos(psi);
w = -V_relative*sin(gamma);

%To speed up calculations, only perform sine and cosine calcs once
sin_gamma = sin(gamma);
cos_gamma = cos(gamma);
tan_gamma = sin_gamma/cos_gamma;
sin_psi = sin(psi);
cos_psi = cos(psi);
sin_lat = sin(lat);
cos_lat = cos(lat);
tan_lat = sin_lat/cos_lat;

%% Sensed acceleration and gyro rates (a(1:3) = accel; a(4:6) = gyro rate)
try
a = interp1q(inputs.time_gyro,inputs.gyro,t);    
% a = specialInterp(t,inputs.gyro,inputs.time_gyro);
catch
    keyboard
end
gyro = a(:);
% Gravitation acceleration; inputs.grav = grav model flag
% grav = gravity(x(1),x(2),x(3),inputs.mu,inputs.J2,inputs.grav);
% g_local_r = norm(grav);
g_local_r = inputs.mu/r^2;

% Determine equivalent Euler angles (in rad)
[~,~,roll] = Quat2Eul(q0,q1,q2,q3);
bank_angle = roll;

% Planet rotation rate
omega_planet = inputs.omegaplanet;

%% Aerodynamic coefficients
% Local speed of sound (m/s)
speedSound = interp1q(inputs.actual.time,inputs.actual.cs,t);

% Calculate Mach infinity (or Mach vehicle when W = 0)
mach = V_relative/speedSound;

% Rotation vector from veh. to body
Rveh2body = R_veh2body_wq(q0,q1,q2,q3);

% Vehicle velocity wrt ground (and wind when W = 0) in body frame
vel_body = Rveh2body*[u;v;w];

% Calculate angle-of-attack (deg) and sideslip angle (deg)
[alpha,beta] = Vw_2_alphabeta(vel_body);
alpha = deg2rad(alpha);
beta = deg2rad(beta);
[CL,CD] = aero_interp2_CLCD(inputs.AOA_grid,inputs.Mach_grid,inputs.data_grid,alpha,beta,mach);

A = inputs.RefArea;

if length(inputs.mass) == 1
    m = inputs.mass;
else
    m = interp1q(inputs.CDtime,inputs.mass,t);
end
% m = inputs.mass;

% Calculate acceleration due to lift and drag
drag = (1/2)*rho*(V_relative^2)*CD*A;
accel_drag = drag/m;

lift = (1/2)*rho*(V_relative^2)*CL*A;
accel_lift = lift/m;

%% Equations of motion
r_dot          = V_relative*sin_gamma;
lat_dot   = V_relative*cos_gamma*sin_psi/r ;
longitude_dot  = V_relative*cos_gamma*cos_psi / (r*cos_lat);
V_relative_dot = -accel_drag - g_local_r*sin_gamma + omega_planet^2*r*cos_lat*(sin_gamma*cos_lat-cos_gamma*sin_lat*sin_psi);
gamma_dot      = (1/V_relative) * (accel_lift*cos(bank_angle) - g_local_r*cos_gamma + (V_relative)^2*cos_gamma/r + ...
    2*omega_planet*V_relative*cos_lat*cos_psi + (omega_planet^2)*r*cos_lat*(cos_gamma*cos_lat+sin_gamma*sin_lat*sin_psi));
psi_dot        = (1/V_relative) * ( accel_lift*sin(bank_angle)/cos_gamma - ((V_relative^2)/r)*cos_gamma*cos_psi*tan_lat + ...
    2*omega_planet*V_relative*(tan_gamma*cos_lat*sin_psi-sin_lat) - ((omega_planet^2)*r/cos_gamma)*sin_lat*cos_lat*cos_psi );
effective_omega = gyro-(1/r)*Rveh2body*[v;-u;-v*tan(lat)];
qdot = 0.5.*qmult_left(q0,q1,q2,q3)*[0;effective_omega];
Pinfdot = rho*g_local_r*w;
rhoinfdot = (rho^2)*g_local_r*w/Pinf;
xdot = [r_dot;lat_dot;longitude_dot;V_relative_dot;gamma_dot;psi_dot;qdot;Pinfdot;rhoinfdot];

if sum(isnan(xdot))>=1
    keyboard
end

return
