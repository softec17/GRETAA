function F = SixDOFEntryEOMVinhEOMjacMSLEOMAccel(t,x,inputs)
% Jacobian of the equations of motion with respect to the state variables

% Length of states
num_x = length(x);
% Preallocate space for the jacobian
F = zeros(num_x);

%% Define state variables
r          = x(1);              % m,   planetocentric radius
lat   = x(2);                   % rad, geocentric lat
lon = x(3);
u = x(4);   
v      = x(5); 
w        = x(6);  
q_J2000_2_DS_0 = x(7);          % Quaternions
q_J2000_2_DS_1 = x(8);
q_J2000_2_DS_2 = x(9);
q_J2000_2_DS_3 = x(10);
q_J2000_2_MCMF_0 = x(11);          
q_J2000_2_MCMF_1 = x(12);
q_J2000_2_MCMF_2 = x(13);
q_J2000_2_MCMF_3 = x(14);
rho = x(15);
Pinf = x(16);

% Planet rotation rate
omega_planet = inputs.omegaplanet;

% Gravity terms
mu = inputs.mu;
J2 = inputs.J2;

%% Rotation matrices
q_J2000_DS = [q_J2000_2_DS_0,q_J2000_2_DS_1,q_J2000_2_DS_2,q_J2000_2_DS_3];
q_J2000_MCMF = [q_J2000_2_MCMF_0,q_J2000_2_MCMF_1,q_J2000_2_MCMF_2,q_J2000_2_MCMF_3];
R_J_DS = quat2dcm(q_J2000_DS); R_DS_2_J = R_J_DS';
R_J_MCMF = quat2dcm(q_J2000_MCMF);
R_MCMF_LH = R2(3*pi/2-lat)*R3(lon);
R_DS_LH = R_MCMF_LH*R_J_MCMF*R_DS_2_J; R_LH_DS = R_DS_LH';

%% Acceleration
% Acceleration
acc_list = inputs.accel;
acc_time_list = inputs.time_accel;
acc = interp1q(acc_time_list,acc_list,t); acc = acc(:);

% IMU converted to LH frame
acc = R_DS_LH*acc;
ax = acc(1);
ay = acc(2);
az = acc(3);

% Angular acceleration
gyro = interp1q(inputs.time_gyro,inputs.gyro,t); gyro = gyro(:);    % In DS frame
om_x = gyro(1); om_y = gyro(2); om_z = gyro(3);

%% Position and velocity equation jacobians
%drdot
F(1,6) = -1;        
%dlatdot
F(2,1) = -u/r^2;    
F(2,4) = 1/r;       
%dlongdot
F(3,1) = -v/(r^2*cos(lat)); 
F(3,2) = v*tan(lat)/(r*cos(lat));
F(3,4) = 1/(r*cos(lat));
%dudot
F(4,1) = -(1/r^2)*(u*w-v^2*tan(lat));
F(4,2) = -v^2/(r*(cos(lat))^2);
F(4,4) = w/r;
F(4,5) = -2*v*tan(lat)/r;
F(4,6) = u/r;
%dvdot
F(5,1) = -(1/r^2)*(u*v*tan(lat)+v*w);
F(5,2) = u*v/(r*(cos(lat))^2);
F(5,4) = v*tan(lat)/r;
F(5,5) = (u*tan(lat)+w)/r;
F(5,6) = v/r;
%dwdot
F(6,1) = (u^2+v^2)/r^2;
F(6,4) = -2*u/r;
F(6,5) = -2*v/r;

%% Sensed acceleration terms

xx = [lat;lon;q_J2000_2_DS_0;q_J2000_2_DS_1;q_J2000_2_DS_2;q_J2000_2_DS_3;q_J2000_2_MCMF_0;q_J2000_2_MCMF_1;q_J2000_2_MCMF_2;q_J2000_2_MCMF_3];
delh = 1e-11;
dV_dx = zeros(3,10);
for ii = 1:length(xx)
    xplus = xx; xminus = xx;
    dx = abs(xx(ii))*delh + eps;
    xplus(ii) = xplus(ii) + dx;
    xminus(ii) = xminus(ii) - dx;
    xdotplus = AccTimeRate(xplus,acc);
    xdotminus = AccTimeRate(xminus,acc);
    dV_dx(:,ii) = (xdotplus-xdotminus)/(2*dx);
end
F(4:6,:) = F(4:6,:) + [zeros(3,1),dV_dx(:,1:2),zeros(3,3),dV_dx(:,3:10),zeros(3,2)];

%% Gravity terms
dg_dx = [[                        (6*J2*mu*sin(2*lat))/r^5,        -(3*J2*mu*cos(2*lat))/r^4, 0]
[                                               0,                                0, 0]
[ - (2*mu)/r^3 - (6*J2*mu*(3*cos(lat)^2 - 2))/r^5, -(9*J2*mu*cos(lat)*sin(lat))/r^4, 0]];
F(4:6,1:3) = F(4:6,1:3) + dg_dx;

%% Quaternion states

dq_J2000_2_DS = [[      0, -om_x/2, -om_y/2, -om_z/2]
[ om_x/2,       0,  om_z/2, -om_y/2]
[ om_y/2, -om_z/2,       0,  om_x/2]
[ om_z/2,  om_y/2, -om_x/2,       0]];
F(7:10,7:10) = dq_J2000_2_DS;

dq_J2000_2_MCMF = [[              0,               0,              0, -omega_planet/2]
[              0,               0, omega_planet/2,               0]
[              0, -omega_planet/2,              0,               0]
[ omega_planet/2,               0,              0,               0]];
F(11:14,11:14) = dq_J2000_2_MCMF;

%% Atmosphere terms
drho_dx = [ -(2*mu*rho^2*w)/(Pinf*r^3), 0, 0, 0, 0, (mu*rho^2)/(Pinf*r^2), 0, 0, 0, 0, 0, 0, 0, 0, (2*mu*rho*w)/(Pinf*r^2), -(mu*rho^2*w)/(Pinf^2*r^2)];
dPinf_dx = [ -(2*mu*rho*w)/r^3, 0, 0, 0, 0, (mu*rho)/r^2, 0, 0, 0, 0, 0, 0, 0, 0, (mu*w)/r^2, 0];
F(15:16,:) = [drho_dx;dPinf_dx];

%% Numerically evaluate jacobian
% Fnumeric = zeros(num_x);
% delh = 1e-11;
% parfor jj = 1:num_x
%     xplus = x; xminus = x;
%     dx = abs(x(jj))*delh + eps;
%     xplus(jj) = xplus(jj) + dx;
%     xminus(jj) = xminus(jj) - dx;
%     xdotplus = SixDOFEntryEOMVinhEOMMSLEOMAccel(t,xplus,inputs);
%     xdotminus = SixDOFEntryEOMVinhEOMMSLEOMAccel(t,xminus,inputs);
%     Fnumeric(:,jj) = (xdotplus-xdotminus)/(2*dx);
% end
% %
% keyboard

% F = Fnumeric;
%%
% if t>150
%    keyboard
% end

return