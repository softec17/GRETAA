function [alpha,beta] = vec_Vw_2_alphabeta(Vwvec)
% Takes a vector of veh. velocity wrt to wind in body frame and finds
% angle-of-attack (alpha) and sideslip angle (beta)
% Can handle a matrix with various velocity vectors as submatrices 
%
% Input
% Vwvec = velocity column vector of vehicle vel. wrt to the wind in body-fixed
%      frame
%
% Outputs
% alpha = angle-of-attack, deg (also alpha_x as opposed to total AOA)
% beta = sideslip angle, deg
%
% Assumption: Vwvec is a column matrix

% Length of column matrix of vectors
N = length(Vwvec);
% Preallocate space for alpha and beta
alpha = zeros(N,1);
beta = zeros(N,1);
% Calculate alpha and beta for each velocity vector
for ii = 1:1:N
   [alphahold,betahold] = Vw_2_alphabeta(Vwvec(ii,:));
   alpha(ii) = alphahold;
   beta(ii) = betahold;
end

return