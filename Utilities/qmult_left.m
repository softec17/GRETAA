function qmat = qmult_left(q0,q1,q2,q3)
% Returns the left matrix needed for quaternion multiplication
% Inputs
% q0,q1,q2,q3 = quaternions (q0 = scalar)
% Output
% qmat = left matrix

qmat = [q0 -q1 -q2 -q3;
        q1  q0 -q3  q2;
        q2  q3  q0 -q1;
        q3 -q2  q1  q0;];

return