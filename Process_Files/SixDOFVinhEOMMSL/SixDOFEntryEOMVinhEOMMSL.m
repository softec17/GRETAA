function xdot = SixDOFEntryEOMVinhEOMMSL(t,x,inputs)
% Equations of motions for 6-DOF where position is r, lat, lon and
% velocity is v, gamma(FPA) and psi (heading angle) from the local
% horizontal frame

disp(t)

%% Define state variables
r          = x(1);  %m,   planetocentric radius
lat   = x(2);  %rad, geocentric lat
lon = x(3);    %rad, longitude     
V_relative = x(4);  %m/s, relative velocity
gamma      = x(5);  %rad, relative flight path angle
psi        = x(6);  %rad, relative heading angle
qJDS0 = x(7);          % Quaternions J2000_DS
qJDS1 = x(8);
qJDS2 = x(9);
qJDS3 = x(10);
% Pinf = x(11);
% rho = x(12);
% qJMCMF0 = x(13);
% qJMCMF1 = x(14);
% qJMCMF2 = x(15);
% qJMCMF3 = x(16);
% rho = interp1q(inputs.rhoTime,inputs.rhoHist,t);
alt = r - inputs.Re;

qJMCMF0 = x(11);
qJMCMF1 = x(12);
qJMCMF2 = x(13);
qJMCMF3 = x(14);

qJDS = [qJDS0,qJDS1,qJDS2,qJDS3];
qJMCMF = [qJMCMF0,qJMCMF1,qJMCMF2,qJMCMF3];
R_J_DS = quat2dcm(qJDS); R_DS_2_J = R_J_DS';
R_J_MCMF = quat2dcm(qJMCMF);
R_MCMF_LH = R2(3*pi/2-lat)*R3(lon);
R_DS_LH = R_MCMF_LH*R_J_MCMF*R_DS_2_J; R_LH_DS = R_DS_LH';
qLHDS = dcm2quat(R_LH_DS);

% Other secondary states
u = V_relative*cos(gamma)*sin(psi);
v = V_relative*cos(gamma)*cos(psi);
w = -V_relative*sin(gamma);

% To speed up calculations, only perform sine and cosine calcs once
sin_gamma = sin(gamma);
cos_gamma = cos(gamma);
tan_gamma = sin_gamma/cos_gamma;
sin_psi = sin(psi);
cos_psi = cos(psi);
sin_lat = sin(lat);
cos_lat = cos(lat);
tan_lat = sin_lat/cos_lat;

%% Sensed acceleration and gyro rates (a(1:3) = gyro rate)
try
    a = interp1q(inputs.time_gyro,inputs.gyro,t);    
catch
    keyboard
end
gyro = a(:);    % In DS frame

% Gravitation acceleration; inputs.grav = grav model flag
grav = gravity(x(1),x(2),x(3),inputs.mu,inputs.J2,inputs.grav);
g_local_r = norm(grav);
% g_local_r = inputs.mu/r^2;

% Determine equivalent Euler angles (in rad)
[~,~,roll] = Quat2Eul(qLHDS(1),qLHDS(2),qLHDS(3),qLHDS(4));
bank_angle = roll;

% Planet rotation rate
omega_planet = inputs.omegaplanet;

%% Aerodynamic coefficients
% Local speed of sound (m/s)
if alt <= 127878.107345943 && alt >= 2781.16091622789
    speedSound = interp1(inputs.PressureInputs.speedSound.altlist,inputs.PressureInputs.speedSound.speedsoundlist,alt);
else
    if alt > 127878.107345943 
        speedSound = inputs.PressureInputs.speedSound.speedsoundlist(1);
    else
        speedSound = inputs.PressureInputs.speedSound.speedsoundlist(end);
    end
end

% Calculate Mach infinity (or Mach vehicle when W = 0)
mach = V_relative/speedSound;

% Rotation vector from veh. to body
Rveh2body = R_LH_DS;

% Vehicle velocity wrt ground (and wind when W = 0) in body frame
vel_body = Rveh2body*[u;v;w];

% Calculate angle-of-attack (deg) and sideslip angle (deg)
[alpha,beta] = Vw_2_alphabeta(vel_body);
alpha = deg2rad(alpha);
beta = deg2rad(beta);
try
[CL,CD] = aero_interpMSL2_CLCD(inputs.AOA_grid,inputs.Mach_grid,inputs.data_grid,alpha,beta,mach);
catch
    keyboard
end

A = inputs.RefArea;

if length(inputs.mass) == 1
    m = inputs.mass;
else
    m = interp1q(inputs.mass_time,inputs.mass,t);
end

% Calculate acceleration due to lift and drag
rho_use = interp1(inputs.rhoAlt,inputs.rhoHist,alt);

drag = (1/2)*rho_use*(V_relative^2)*CD*A;
accel_drag = drag/m;

lift = (1/2)*rho_use*(V_relative^2)*CL*A;
accel_lift = lift/m;

%% Equations of motion
r_dot     = V_relative*sin_gamma;
lat_dot   = V_relative*cos_gamma*sin_psi/r ;
longitude_dot  = V_relative*cos_gamma*cos_psi / (r*cos_lat);
V_relative_dot = -accel_drag - g_local_r*sin_gamma + omega_planet^2*r*cos_lat*(sin_gamma*cos_lat-cos_gamma*sin_lat*sin_psi);
gamma_dot      = (1/V_relative) * (accel_lift*cos(bank_angle) - g_local_r*cos_gamma + (V_relative)^2*cos_gamma/r + ...
    2*omega_planet*V_relative*cos_lat*cos_psi + (omega_planet^2)*r*cos_lat*(cos_gamma*cos_lat+sin_gamma*sin_lat*sin_psi));
psi_dot        = (1/V_relative) * ( accel_lift*sin(bank_angle)/cos_gamma - ((V_relative^2)/r)*cos_gamma*cos_psi*tan_lat + ...
    2*omega_planet*V_relative*(tan_gamma*cos_lat*sin_psi-sin_lat) - ((omega_planet^2)*r/cos_gamma)*sin_lat*cos_lat*cos_psi );
% Pinfdot = rho*g_local_r*w;
% rhoinfdot = (rho^2)*g_local_r*w/Pinf;

q_J2000_2_DS_dot = quatEOM_general(t,qJDS,gyro);
q_J2000_2_MCMF_dot = quatEOM(t,qJMCMF,omega_planet);

% xdot = [r_dot;lat_dot;longitude_dot;V_relative_dot;gamma_dot;psi_dot;q_J2000_2_DS_dot;Pinfdot;rhoinfdot;q_J2000_2_MCMF_dot];
xdot = [r_dot;lat_dot;longitude_dot;V_relative_dot;gamma_dot;psi_dot;q_J2000_2_DS_dot;q_J2000_2_MCMF_dot];

if sum(isnan(xdot))>=1
    keyboard
end
% keyboard
return
