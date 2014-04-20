function [x2estimate, X2estimate, Uncertainty, P2estimate,zk,Kk,Qk,...
    qk_list,Deltaq_k_list,qk_mean_list,Qk_list] = AdaptiveFilter2(t2,X1,P1,t1,...
    measurements,HMatrixFcnHandle,RMatrixFcnHandle,...
    inputs,error_level,kk,N,qk_list,Deltaq_k_list,qk_mean_list,Qk_list,Qk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION: AdaptiveFilter
%BY: Soumyo Dutta
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
%   zk - innovation matrix
%   Kk - Kalman Gain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% State propogation
len_state = length(X1);         % Size of covariance
% Propagating state and covariance
if kk > N       % Adaptive filter choice for process noise covariance
    Qkuse = Qk;                     % covEOM will use this Q instead of Q from QMat function
    kala = 1;
else
    Qkuse = zeros(len_state);       % Send in zero matrix -- covEOM will calculate Q
    kala = 2;
end
[X2nominal,P2nominal,P2nominal_noQ] = propagateStateAdaptiveFilter(X1,P1,t1,t2,inputs,kk,N,Qkuse);
delX = 1;
iter = 1;
iter_max = 1;
while delX>inputs.reltol&&iter<=iter_max
    %% Kalman update (i.e. measurement update)
    % Calculate measurement noise matrix
    Rk = RMatrixFcnHandle(t2,X2nominal,measurements,inputs,error_level);
    try
        [Hk,zk] = HMatrixFcnHandle(t2,X2nominal,measurements,inputs);
    catch
        keyboard;
    end
    A = P2nominal*Hk';
    % keyboard
    B = Hk*P2nominal*Hk'+Rk;
    Kk = A/B;
    cond_num = rcond(B);
    % fprintf('Current time: %8.4f sec\n',t2)     % Prints current observation time to screen
    disp(['Current time: ',num2str(t2),' sec. Condition number: ',num2str(cond_num),'. ',num2str(kala),' Iter.: ',num2str(iter)]);
%     disp(['Current time: ',num2str(t2),' sec. Condition number: ',num2str(cond_num),'. ',num2str(kala)]);

    % Covariance update
    Resid = (eye(len_state)-Kk*Hk);
    P2estimate = Resid*P2nominal*Resid'+Kk*Rk*Kk';
    % Deviation in state estimate
    x2estimate = Kk*zk;
    % State vector update
    X2estimate = X2nominal + x2estimate;
    % Iterative EKF tolerance calculations and resetting states if
    % iteration needed
    delX = max(abs(x2estimate));
%     keyboard;
    X2nominal = X2estimate;
    P2nominal = P2estimate;
    iter = iter + 1;
%     keyboard
    
%     act_rad = interp1q(inputs.actual.time,inputs.actual.gcrad,t2);
%     if abs(act_rad-X2estimate(1))./1e3 > 1e3
%         keyboard
%     end
end
%Calculate error in each element
Uncertainty = sqrt(diag(P2estimate))';

%% Adaptive filter covariance matching
% Calculate process noise vector approx.
qk = X2estimate - X2nominal;
Deltaq_k = P2nominal_noQ - P2estimate;

% Store calculated values for future use
qk_list{kk} = qk;
Deltaq_k_list{kk} = Deltaq_k;

% Calculate process noise for next step -- using Myers/Tapley cumulation
if kk > N
    qk_mean = qk_mean_list{kk-1};
    qk_mean = qk_mean + (1/N)*(qk - qk_list{kk-N});
    Qk = Qk + (1/(N-1))*((qk-qk_mean)*(qk-qk_mean)'-(qk_list{kk-N}-qk_mean)*(qk_list{kk-N}-qk_mean)'+...
        + (1/N)*(qk - qk_list{kk-N})*(qk - qk_list{kk-N})' + ((N-1)/N)*(Deltaq_k_list{kk-N}-Deltaq_k));
else        % Sample size for q not reached yet -- no change in Q, store qk_mean
    if kk ~= 1
        if kk ~= N
            qk_mean = qk_mean_list{kk-1};
            qk_mean = (qk_mean*(kk-1) + qk)/kk;
            Qk = zeros(len_state);
        else
            qk_mean = qk_mean_list{kk-1};
            qk_mean = (qk_mean *(kk-1)+qk)/kk;
            Qk = zeros(len_state);
            for ii = 1:N
                Qk = (qk - qk_mean)*(qk - qk_mean)' - ((N-1)/N)*Deltaq_k;
            end
            Qk = (1/(N-1)).*Qk;
        end
    else
        qk_mean = qk;
        Qk = zeros(len_state);
    end
end

% % Calculate process noise for next step -- using actual Myers/Tapley
% % equations not the cumulation equations
% if kk > N
%     qk_mean = zeros(size(qk));
%     Qk = zeros(
%     for jj = 1:N
%         qk_mean = qk_mean_list{kk-1};
%         qk_mean = qk_mean + (1/N)*(qk - qk_list{kk-N});
%         Qk = Qk + (1/(N-1))*((qk-qk_mean)*(qk-qk_mean)'-(qk_list{kk-N}-qk_mean)*(qk_list{kk-N}-qk_mean)'+...
%             + (1/N)*(qk - qk_list{kk-N})*(qk - qk_list{kk-N})' + ((N-1)/N)*(Deltaq_k_list{kk-N}-Deltaq_k));
%     end
% else
%
% end

% Check for positive definiteness of the Q-matrix
[~,eigenval] = eig(Qk);
eigenval_vec = diag(eigenval);
ind_neg_eig = find(eigenval_vec<=0,1);
if isempty(ind_neg_eig)~=1
    Qk_calc = Qk;
    for uu = 1:length(eigenval_vec)
        Qk(uu,uu) = abs(Qk_calc(uu,uu));
    end
end

% % Make position vector process noise and its covariance zero
% for ww = 1:3
%     Qk(ww,:) = zeros(1,len_state);
%     Qk(:,ww) = zeros(len_state,1);
% end

qk_mean_list{kk} = qk_mean;
Qk_list{kk} = Qk;

return
