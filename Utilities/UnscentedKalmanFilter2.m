function [x2estimate, X2estimate, Uncertainty, P2estimate] = ...
    UnscentedKalmanFilter2(t2,X1,P1,t1,...
    measurements,HMatrixFcnHandle,RMatrixFcnHandle,QMatrixFcnHandle,...
    inputs,error_level,kappa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UnscentedKalmanFilter
%
% Inputs
% t2 = time of the next measurement update
% X1 = best estimate of state at t = t1
% P1 = best estimate of the state covariance at t = t1
% t1 = time of the previous best estimate
% measurements = vector of measurements at t = t2
% HMatrixFcnHandle = handle for the meas. jacobian matix (jacobian is not
% used for UKF but the function also calculates the predicted meas. which 
% is used in the calc.
% RMatrixFcnHandle = measurement uncertainty matrix
% QMatrixFcnHandle = process noise uncertainty matrix
%
% Outputs
% x2estimate = state deviation in the filter
% X2estimate = updated state estimate (after measurement update)
% Uncertainty = 1 sigma uncertanty of the states (no correlations)
% P2estimate = update covariance matrix (after measurement update)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prints current observation time to screen
fprintf('Current time: %8.4f sec\n',t2)
% State vector size
N = length(X1);

% Measurement vector size
M = length(measurements);

%% Propagate state and covariance to the next point
Beforeprop_xsig = CreateSigmaPt(X1,P1,kappa);
% Beforeprop_xsig = real(Beforeprop_xsig);
if isreal(Beforeprop_xsig)~=1
    keyboard;
end
Afterprop_xsig = zeros(size(Beforeprop_xsig));
EOM = inputs.EOM;
options = odeset('RelTol',1e-11);
% Propagate signma points
parfor ii = 1:2*N+1
    [t,x] = ode45((@(t,x) EOM(t,x,inputs)),[t1 t2],Beforeprop_xsig(:,ii),options);
    Afterprop_xsig(:,ii) = x(end,:)';
end
if isreal(Afterprop_xsig)~=1
    keyboard;
end
[X2nominal,P2nominal] = UndoSigmaPt(Afterprop_xsig,kappa);
if isreal(X2nominal)~=1
    keyboard;
end
if isreal(P2nominal)~=1
    keyboard;
end
% X2nominal = real(X2nominal);
% P2nominal = real(P2nominal);
P2nominal = P2nominal + QMatrixFcnHandle(t1,X1,inputs);
if isreal(P2nominal)~=1
    keyboard;
end
%% Kalman measurement update
% We need new sigma points and cannot use the last sigma points after the
% time propagation. That is because the distribution of the sigma points
% are no longer at the preferred weighting points that are desired
Beforemeas_xsig = CreateSigmaPt(X2nominal,P2nominal,kappa);
if isreal(Beforemeas_xsig)~=1
    keyboard;
end
% Beforemeas_xsig = real(Beforemeas_xsig);
Aftermeas_ysig = zeros(M,2*N+1);
yk = zeros(M,1);
Py = zeros(M);
Pxy = zeros(N,M);
alpha = 1e-3;
beta = 2;
lamda = alpha^2*(N+kappa)-N;
W_0_mean = lamda/(N+lamda);
W_0_cov = W_0_mean + 1 - alpha^2 + beta;
W_i = 1/(2*(N+lamda));
% Calculate sigma points after the measurement equation process
parfor ii = 1:2*N+1
    % The measurement equation will also calculate the measurement
    % sensitivity matrix (Hk) but this will not be used
    % zk is the residual between the actual measurement and the prediction
    % We want the predicted value only
    [Hk,zk] = HMatrixFcnHandle(t2,Beforemeas_xsig(:,ii),measurements,inputs);
    zk = -1.*(zk(:) - measurements(:));
    Aftermeas_ysig(:,ii) = zk(:);
    if ii ==1
        yk = yk + W_0_mean.*zk(:);
    else
        yk = yk + W_i.*zk(:);
    end
end
if isreal(Aftermeas_ysig)~=1
    keyboard;
end
% Use the sigma points after going through the measurement process for
% update to the state and covariance
parfor ii = 1:2*N+1
    residy = Aftermeas_ysig(:,ii)-yk;
%     residy = real(residy);
    residx = Beforemeas_xsig(:,ii)-X2nominal;
%     residx = real(residx);
    if ii == 1
        Py = Py + W_0_cov.*(residy*residy');
        Pxy = Pxy + W_0_cov.*(residx*residy');
    else
        Py = Py + W_i.*(residy*residy');
        Pxy = Pxy + W_i.*(residx*residy');
    end
end
Rk = RMatrixFcnHandle(t2,X2nominal,measurements,inputs,error_level);
Py = Py + Rk;

%% Kalman gain calculation and state/covariance update
Kk = Pxy/Py;
x2estimate = Kk*(measurements(:)-yk);
X2estimate = X2nominal + x2estimate;
P2estimate = P2nominal - Kk*Py*Kk';
Uncertainty = sqrt(diag(P2estimate))';

if isreal(x2estimate)~=1
    keyboard;
end

if isreal(P2estimate)~=1
    keyboard;
end


% x2estimate = real(x2estimate);
% X2estimate = real(X2estimate);
% P2estimate = real(P2estimate);
% Uncertainty = real(Uncertainty);

return
