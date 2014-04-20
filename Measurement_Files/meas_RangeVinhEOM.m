function [range,azimuth,elevation] = meas_Range(t,x,inputs)
% Inputs
% t = time (sec)
% x = states
% inputs = structure containing information for state propagation
%
% Output
% range = range from station to s/c (m)
% azimuth = heading angle from station to s/c (deg)
% elevation = elevation angle from station to s/c (deg)

% Current states
r = x(1);
lat = x(2);
lon = x(3);
V = x(4);
gamma      = x(5);  %rad, relative flight path angle
psi        = x(6);  %rad, relative heading angle
% Quaternions
q0 = x(7);
q1 = x(8);
q2 = x(9);
q3 = x(10);
% Atmospheric parameters
Pinf = x(11);
rhoinf = x(12);

% Other secondary states
u = V*cos(gamma)*sin(psi);
v = V*cos(gamma)*cos(psi);
w = -V*sin(gamma);

% Station position (assuming non-moving station)
r_station = inputs.Radar.r_station;
t_station = inputs.Radar.t_station;
rstation_t = interp1(t_station,r_station,t);

% S/c position with radar measurement
r_sc = r.*[cos(lat)*cos(lon);cos(lat)*sin(lon);sin(lat)];
delta_R = r_sc - rstation_t;
range = sqrt(delta_R'*delta_R);
azimuth = atan2(delta_R(2),delta_R(1));
elevation = asin2(delta_R(3),range);

return