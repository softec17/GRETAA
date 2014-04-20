function covdot = SixDOFEntryCOVVinhEOMMSLEOMAccel(t,xcov,xstates,inputs)
% Calculates the derivative of the covariance matrix needed for its
% propagation in time.
% t = current time
% xcov = upper triangle covariance matrix elements, row by row
% xstates = states being estimated
% inputs = constants needed for propagation of EOMs

% Jacobian of the equations of motion with respect to the state variables
F = SixDOFEntryEOMVinhEOMjacMSLEOMAccel(t,xstates,inputs);

% Process noise vector associated with the equations of motion
Q = SixDOFEntryEOMQMatVinhEOMMSLEOMAccel(t,xstates,inputs);

% Make the upper triangle covariance vector into the symmetric matrix
len_cov = length(xcov);
xcov1 = xcov(1:len_cov);
cov1 = makemat(xcov1);

% Calculate the covariance rates
covdotmat1 = F*cov1 + cov1'*F'+ Q;

% Make the covariance rate matrix into a vector
covdot = makevec(covdotmat1);

return