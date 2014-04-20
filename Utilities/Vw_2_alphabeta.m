function [alpha,beta] = Vw_2_alphabeta(Vw)
% Calculate angle-of-attack (alpha) and sideslip angle (beta) when given
% the velocity of the vehicle wrt to the wind in body-fixed frame
%
% Input
% Vw = velocity column vector of vehicle vel. wrt to the wind in body-fixed
%      frame
%
% Outputs
% alpha = angle-of-attack, deg (also alpha_x as opposed to total AOA)
% beta = sideslip angle, deg

% Calculate alpha and beta
Vwmag = norm(Vw);
alpha = atan(Vw(3)/Vw(1));
% % beta = atan(Vw(2)/sqrt(Vw(1)^2+Vw(3)^2));
beta = asin(Vw(2)/Vwmag);
% Convert angles to degrees
alpha = rad2deg(alpha);
beta = rad2deg(beta);

return