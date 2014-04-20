function [Ubank,Upitch,Uazimuth] = Quat2EulUncer2(q0list,q1list,q2list,q3list,bank,pitch,azimuth,uq0list,uq1list,uq2list,uq3list)
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

Ubank = zeros(size(q0list));
Upitch = Ubank;
Uazimuth = Ubank;

for ii = 1:length(q0list)
    
    q0 = q0list(ii);
    q1 = q1list(ii);
    q2 = q2list(ii);
    q3 = q3list(ii);
        
    uq0 = uq0list(ii);
    uq1 = uq1list(ii);
    uq2 = uq2list(ii);
    uq3 = uq3list(ii);
    
j_roll = [ 0, 0, 0, 0, 0, 0, ((2*q1)/(q0^2 - q1^2 - q2^2 + q3^2) - (2*q0*(2*q0*q1 + 2*q2*q3))/(q0^2 - q1^2 - q2^2 + q3^2)^2)/((2*q0*q1 + 2*q2*q3)^2/(q0^2 - q1^2 - q2^2 + q3^2)^2 + 1), ((2*q0)/(q0^2 - q1^2 - q2^2 + q3^2) + (2*q1*(2*q0*q1 + 2*q2*q3))/(q0^2 - q1^2 - q2^2 + q3^2)^2)/((2*q0*q1 + 2*q2*q3)^2/(q0^2 - q1^2 - q2^2 + q3^2)^2 + 1), ((2*q3)/(q0^2 - q1^2 - q2^2 + q3^2) + (2*q2*(2*q0*q1 + 2*q2*q3))/(q0^2 - q1^2 - q2^2 + q3^2)^2)/((2*q0*q1 + 2*q2*q3)^2/(q0^2 - q1^2 - q2^2 + q3^2)^2 + 1), ((2*q2)/(q0^2 - q1^2 - q2^2 + q3^2) - (2*q3*(2*q0*q1 + 2*q2*q3))/(q0^2 - q1^2 - q2^2 + q3^2)^2)/((2*q0*q1 + 2*q2*q3)^2/(q0^2 - q1^2 - q2^2 + q3^2)^2 + 1)];
j_pitch = [ 0, 0, 0, 0, 0, 0, (2*q2)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2), -(2*q3)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2), (2*q0)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2), -(2*q1)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2)];
j_yaw = [ 0, 0, 0, 0, 0, 0, ((2*q3)/(q0^2 + q1^2 - q2^2 - q3^2) - (2*q0*(2*q0*q3 + 2*q1*q2))/(q0^2 + q1^2 - q2^2 - q3^2)^2)/((2*q0*q3 + 2*q1*q2)^2/(q0^2 + q1^2 - q2^2 - q3^2)^2 + 1), ((2*q2)/(q0^2 + q1^2 - q2^2 - q3^2) - (2*q1*(2*q0*q3 + 2*q1*q2))/(q0^2 + q1^2 - q2^2 - q3^2)^2)/((2*q0*q3 + 2*q1*q2)^2/(q0^2 + q1^2 - q2^2 - q3^2)^2 + 1), ((2*q1)/(q0^2 + q1^2 - q2^2 - q3^2) + (2*q2*(2*q0*q3 + 2*q1*q2))/(q0^2 + q1^2 - q2^2 - q3^2)^2)/((2*q0*q3 + 2*q1*q2)^2/(q0^2 + q1^2 - q2^2 - q3^2)^2 + 1), ((2*q0)/(q0^2 + q1^2 - q2^2 - q3^2) + (2*q3*(2*q0*q3 + 2*q1*q2))/(q0^2 + q1^2 - q2^2 - q3^2)^2)/((2*q0*q3 + 2*q1*q2)^2/(q0^2 + q1^2 - q2^2 - q3^2)^2 + 1)];
uncer_xvec = [0; 0; 0; 0; 0; 0; uq0; uq1; uq2; uq3];

Ubank(ii) = j_roll*uncer_xvec;
Upitch(ii) = j_pitch*uncer_xvec;
Uazimuth(ii) = j_yaw*uncer_xvec;

end

return