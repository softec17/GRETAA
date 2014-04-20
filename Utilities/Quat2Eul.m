function [azimuth,pitch,bank] = Quat2Eul(q0,q1,q2,q3)
% Takes quaternions and calculates the equivalent Euler angles
% Inputs
% q0,q1,q2,q3 = quaternions with q0 being the "scalar" value
%
% Outputs
% bank,pitch,azimuth = Euler angles (rad) in 1,2,and 3 axis

m11 = 2.*(q1.*q2 + q0.*q3);
m12 = q0.^2 + q1.^2 - q2.^2 - q3.^2;
m21 = -2.*(q1.*q3 - q0.*q2);
m31 = 2.*(q2.*q3 + q0.*q1);
m32 = q0.^2 - q1.^2 - q2.^2 + q3.^2;

bank = atan2(m31,m32);
pitch = asin(m21);
azimuth = atan2(m11,m12);

return

