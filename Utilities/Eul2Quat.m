function [q0,q1,q2,q3] = Eul2Quat(bank,pitch,azimuth)
% Takes Euler angles and calculates the equivalent quaternion
% Inputs
% bank,pitch,azimuth = Euler angles (rad) in 1,2,and 3 axis
%
% Outputs
% q0,q1,q2,q3 = quaternions with q0 being the "scalar" value

theta1 = azimuth;
theta2 = pitch;
theta3 = bank;

t1 = theta1/2;
t2 = theta2/2;
t3 = theta3/2;

q0 = cos(t1)*cos(t2)*cos(t3)+sin(t1)*sin(t2)*sin(t3);
q1 = cos(t1)*cos(t2)*sin(t3)-sin(t1)*sin(t2)*cos(t3);
q2 = cos(t1)*sin(t2)*cos(t3)+sin(t1)*cos(t2)*sin(t3);
q3 = sin(t1)*cos(t2)*cos(t3)-cos(t1)*sin(t2)*sin(t3);

return
