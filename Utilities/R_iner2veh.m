function rotVI = R_iner2veh(lat,lon)
% Rotation vector to take a vector from inertial frame to vehicle-carried
% frame
%
% Input
% lat = latitude of the vehicle, deg (assumed to be geocentric latitude)
% lon = longitude of the vehicle, deg
%
% Output
% rotVI = rotation matrix to go from inertial to vehicle-carried frame

% Convert angles to radians
lat = deg2rad(lat);
lon = deg2rad(lon);
% Prealocate space
rotVI = zeros(3);
% Define rotation matrix
rotVI(1,1) = -sin(lat).*cos(lon);
rotVI(1,2) = -sin(lat).*sin(lon);
rotVI(1,3) = cos(lat);
rotVI(2,1) = -sin(lon);
rotVI(2,2) = cos(lon);
rotVI(3,1) = -cos(lat).*cos(lon);
rotVI(3,2) = -cos(lat).*sin(lon);
rotVI(3,3) = -sin(lat);

return