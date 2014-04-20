function [x2,P2,P2_noQ] = propagateStateAdaptiveFilter(x1,P1,t1,t2,inputs,kk,N,Qk)
% By: Soumyo Dutta
% Uses EOMs to propagate states from the last time to the next time step

% Converts the upper-triangle of the symmetric covariance to a vector
covvec = makevec(P1);
statelength = length(x1);
covlength = length(covvec);


% x1mod = [x1;covvec;covvec]; % New "state" vector that includes the state and covariance elements

% keyboard
% Type of propagation equation
% EOM = inputs.dynamics;

% Propagating state
% options = odeset('RelTol',inputs.reltol);
% [t,x] = ode45((@(t,x) EOM(t,x,inputs,kk,N,Qk)),[t1 t2],x1mod,options);

% % x2 = x(end,1:statelength)'; % Vector of the propagated state variables
% P2vec = x(end,statelength+1:statelength+covlength)'; % Vector of the propagated covariance
% P2 = makemat(P2vec); % Creates the upper-triangle of the eventually symmetric, covariance matrix
% P2vecnoQ = x(end,statelength+covlength+1:end)'; % Vector of the propagated covariance
% P2_noQ = makemat(P2vecnoQ); % Creates the upper-triangle of the eventually symmetric, covariance matrix

options = odeset('RelTol',inputs.reltol);
EOM = inputs.EOM;
[t,x] = ode45((@(t,x) EOM(t,x,inputs)),[t1 t2],x1,options);
x2 = x(end,:)'; % Vector of the propagated state variables

EOM = inputs.dynamics;
x1mod = [x1;covvec;covvec]; % New "state" vector that includes the state and covariance elements
[t,x] = ode45((@(t,x) EOM(t,x,inputs,kk,N,Qk)),[t1 t2],x1mod,options);
P2vec = x(end,statelength+1:statelength+covlength)'; % Vector of the propagated covariance
P2 = makemat(P2vec); % Creates the upper-triangle of the eventually symmetric, covariance matrix
P2vecnoQ = x(end,statelength+covlength+1:end)'; % Vector of the propagated covariance
P2_noQ = makemat(P2vecnoQ); % Creates the upper-triangle of the eventually symmetric, covariance matrix

return