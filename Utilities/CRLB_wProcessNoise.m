function [Pk1] = CRLB_wProcessNoise(tk1,xk1,tk,xk,Pk,F,Q,H,R,inputs)
% Cramer-Rao Lower Bounds with process noise
% Inputs:
% tk1, xk1      time and state at time k+1
% tk, xk, Pk    time, state, and covariance at time k
% F, Q, H, R    Jacobian handles (state sensitivity, process noise,
%               measurement sensitivity, measurement noise)
% inputs        struct with other propagation information

Fk = F(tk,xk,inputs);
Qk = Q(tk,xk,inputs);
Hk1 = H(tk1,xk1,inputs);
Rk1 = R(tk1,xk1,inputs);
Qinv = inv(Qk);

Dk11 = Fk'*Qinv*Fk;
Dk12 = -Fk'*Qinv;
Dk22 = Qinv + Hk1'/Rk1*Hk1;

Jk1 = Dk22 - Dk12'/(Pk + Dk11)*Dk12;
Pk1 = inv(Jk1);

return