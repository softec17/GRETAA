function rotBV = R_veh2body(bank,pitch,azimuth)
% Rotation matrix to take vectors from the vehicle frame to the body frame
%
% Inputs
% bank = bank or roll angle (rotation angle in x-axis), deg
% pitch = pitch angle (rotation angle in y-axis), deg
% azimuth = azimuth or yaw angle (rotation angle in z-axis), deg
%
% Output
% rotBV = rotation matrix to take matrices from vehicle frame to body frame

% Convert angles from degrees to radians
bank = deg2rad(bank);
pitch = deg2rad(pitch);
azimuth = deg2rad(azimuth);

% Preallocate matrices
rotBV = zeros(3);

rotBV(1,1) = cos(pitch)*cos(azimuth);
rotBV(1,2) = cos(pitch)*sin(azimuth);
rotBV(1,3) = -sin(pitch);
rotBV(2,1) = sin(bank)*sin(pitch)*cos(azimuth) - cos(bank)*sin(azimuth);
rotBV(2,2) = sin(bank)*sin(pitch)*sin(azimuth) + cos(bank)*cos(azimuth);
rotBV(2,3) = sin(bank)*cos(pitch);
rotBV(3,1) = cos(bank)*sin(pitch)*cos(azimuth) + sin(bank)*sin(azimuth);
rotBV(3,2) = cos(bank)*sin(pitch)*sin(azimuth) - sin(bank)*cos(azimuth);
rotBV(3,3) = cos(bank)*cos(pitch);

return