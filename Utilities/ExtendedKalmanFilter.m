function [x2estimate, X2estimate, Uncertainty, P2estimate,zk,Kk] = ...
    ExtendedKalmanFilter(t2,X1,P1,t1,...
    measurements,HMatrixFcnHandle,RMatrixFcnHandle,...
    inputs,error_level)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION: ExtendedKalmanFilter
%BY: Soumyo Dutta and started by John A. Christian 
%LAST MODIFIED: September 2009
%
%This function implements an Extended Kalman Filter (EKF).
%
%INPUTS:
%   t_nominal - vector of times for nominal trajectory (seconds)
%   x_last - state vector for the last time step (mixed units)
%   P_last - vector containting the elements of the last covariance matrix
%   (mixed units)
%   t_last - last time step (seconds)
%   measurements - observations for the current time step
%   HMatrixFcnHandle - function handle that points to H matrix
%   RMatrixFcnHandle - function handle that points to R matrix
%   QMatrixFcnHandle - function handle that points to Q matrix
%   inputs - input structure for preset data (see initialize_states.m)
%
%OUTPUTS:
%   t_history - vector of times for nominal trajectory (seconds)
%   x_n_history - time-history of nominal state vector (mixed units) --
%                 pre-measurement values
%   xk_plus_history - time-history of error in state vector (mixed units)
%   x_updated - time-history of best estimated state vector (mixed units)
%   x_uncertainty - time-history of standard deviation for each state for the
%                   best estimated trajectory (mixed units)
%   P_end - final covariance matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prints current observation time to screen
fprintf('Current time: %8.4f sec\n',t2)
% Size of covariance
len_state = length(X1);
% Propagating state and covariance
[X2nominal,P2nominal] = propagateState(X1,P1,t1,t2,inputs);
% Calculate measurement noise matrix
Rk = RMatrixFcnHandle(t2,X2nominal,measurements,inputs,error_level);
[Hk,zk] = HMatrixFcnHandle(t2,X2nominal,measurements,inputs);
% Kalman gain
% Kk = (P2nominal*Hk')/(Hk*P2nominal*Hk'+Rk);

A = P2nominal*Hk';
B = Hk*P2nominal*Hk'+Rk;
Kk = A/B;
% Kk = (P2nominal*Hk')/(Hk*P2nominal*Hk'+Rk);

% Covariance update
Resid = (eye(len_state)-Kk*Hk);
P2estimate = Resid*P2nominal*Resid'+Kk*Rk*Kk';
% Deviation in state estimate
x2estimate = Kk*zk;     
% State vector update
X2estimate = X2nominal + x2estimate;
%Calculate error in each element
Uncertainty = sqrt(diag(P2estimate))';
% keyboard
return
