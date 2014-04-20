function [roll,pitch,azimuth] = rot321_angles(rotMat)
% Takes a 321 angle rotation matrix and calculates the Euler rotation
% angles from it
%
% Inputs
% rotMat = rotation matrix
%
% Outputs
% Roll = angle is the x axis
% Pitch = angle in the y axis
% Azimuth = angle in the z axis
% Angles in degrees

% [c2c3         -c2s3           s2
% c1s3+c3s1s2   c1c3-s1s2s3     -c2s1
% s1s3-c1c3s2   c1s2s3+c3s1     c1c2]

roll = atan2(rotMat(2,3),rotMat(3,3));
pitch = atan(rotMat(1,3));
azimuth = atan2(rotMat(1,2),rotMat(1,1));

return