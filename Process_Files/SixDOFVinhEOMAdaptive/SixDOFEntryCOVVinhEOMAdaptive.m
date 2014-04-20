function covdot = SixDOFEntryCOVVinhEOMAdaptive(t,xcov,xstates,inputs,kk,N,Qk)
% Calculates the derivative of the covariance matrix needed for its
% propagation in time.
% t = current time
% xcov = upper triangle covariance matrix elements, row by row
% xstates = states being estimated
% inputs = constants needed for propagation of EOMs

% Jacobian of the equations of motion with respect to the state variables
F = SixDOFEntryEOMVinhEOMjacAdaptive(t,xstates,inputs);

% Process noise vector associated with the equations of motion
if kk > N
    Q = Qk;
else
    Q = SixDOFEntryEOMQMatVinhEOMAdaptive(t,xstates,inputs);
end

% Make the upper triangle covariance vector into the symmetric matrix
len_cov = length(xcov)/2;
xcov1 = xcov(1:len_cov);
xcov2 = xcov(len_cov+1:end);
cov1 = makemat(xcov1);
cov2 = makemat(xcov2);

% Calculate the covariance rates
covdotmat1 = F*cov1 + cov1'*F'+ Q;
covdotmat2 = F*cov2 + cov2'*F';

% Make the covariance rate matrix into a vector
covdot = [makevec(covdotmat1);makevec(covdotmat2)];

return