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
u = x(4); % Velocity in North direction (of NED frame)
v = x(5); % Velocity in East direction (of NED frame)
w = x(6); % Velocity in Down direction (of NED frame)
% Quaternions
qJDS0 = x(7);
qJDS1 = x(8);
qJDS2 = x(9);
qJDS3 = x(10);
% Atmospheric parameters
rhoinf = x(15);
Pinf = x(16);


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