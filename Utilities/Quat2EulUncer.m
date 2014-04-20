function [Ubank,Upitch,Uazimuth] = Quat2EulUncer(q0,q1,q2,q3,bank,pitch,azimuth,uq0,uq1,uq2,uq3)
% Calculate uncertainty in the aerospace Euler angles based on the
% uncertainty in the quaternions
%
% Inputs
% q0,q1,q2,q3 = quaternions with q0 being the "scalar" value
% bank,pitch,azimuth = Euler angles (rad) in 1,2,and 3 axis
% uq0,uq1,uq2,uq3 = uncertainty in quaternion
%
% Outputs
% Ubank,Upitch,Uazimuth = uncertainties in Euler angles (rad) in 1,2,and 3
% axis

% For note, the equations based on:
% m11 = 2.*(q1.*q2 + q0.*q3);
% m12 = q0.^2 + q1.^2 - q2.^2 - q3.^2;
% m21 = -2.*(q1.*q3 - q0.*q2);
% m31 = 2.*(q2.*q3 + q0.*q1);
% m32 = q0.^2 - q1.^2 - q2.^2 + q3.^2;
% 
% bank = atan2(m31,m32);
% pitch = asin(m21);
% azimuth = atan2(m11,m12);

m11 = 2.*(q1.*q2 + q0.*q3);
m12 = q0.^2 + q1.^2 - q2.^2 - q3.^2;
m21 = -2.*(q1.*q3 - q0.*q2);
m31 = 2.*(q2.*q3 + q0.*q1);
m32 = q0.^2 - q1.^2 - q2.^2 + q3.^2;

dm11 = 2.*(q2.*uq1+q1.*uq2+q3.*uq0+q0.*uq3);
dm12 = 2.*q0.*uq0 + 2.*q1.*uq1 - 2.*q2.*uq2 - 2.*q3.*uq3;
dm21 = -2.*(q1.*uq3 + q3.*uq1 - q0.*uq2 - q2.*uq0);
dm31 = 2.*(q2.*uq3+q3.*uq2+q0.*uq1+q1.*uq0);
dm32 = 2.*q0.*uq0 - 2.*q1.*uq1 - 2.*q2.*uq2 + 2.*q3.*uq3;

Uazimuth = (1./(sec(azimuth).^2)).*(dm11./m12-(m11./m12.^2).*dm12);
Upitch = (1./cos(pitch)).*(-dm21);
Ubank = (1./(sec(bank).^2)).*(dm31./m32-(m31./m32.^2).*dm32);

Uazimuth = real(Uazimuth);
Upitch = real(Upitch);
Ubank = real(Ubank);

return