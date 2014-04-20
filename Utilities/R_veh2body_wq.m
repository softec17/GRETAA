function rotMat = R_veh2body_wq(q0,q1,q2,q3)
% Rotation matrix from vehicle-fixed frame to body-frame using quaternions
% Inputs
% q0,q1,q2,q3 = quaternions (q0 = scalar)
% Output
% rotMat = rotation matrix from veh. to b-frame

rotMat = [2*q0^2-1+2*q1^2, 2*(q1*q2+q0*q3), 2*(q1*q3-q0*q2);
          2*(q1*q2-q0*q3), 2*q0^2-1+2*q2^2, 2*(q2*q3+q0*q1);
          2*(q1*q3+q0*q2), 2*(q2*q3-q0*q1), 2*q0^2-1+2*q3^2;];

return