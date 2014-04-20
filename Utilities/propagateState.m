function [x2,P2] = propagateState(x1,P1,t1,t2,inputs)
% By: Soumyo Dutta
% Modified: September 2009
% Uses EOMs to propagate states from the last time to the next time step

% Converts the upper-triangle of the symmetric covariance to a vector
covvec = makevec(P1);
covlength = length(covvec);
% New "state" vector that includes the state and covariance elements
x1mod = [x1;covvec];

% Type of propagation equation
EOM = inputs.dynamics;
% keyboard
% Propagating state
% if t1 > t2                          % Backward run -- stiff problem
    options = odeset('RelTol',inputs.reltol);
    % options = odeset('RelTol',inputs.reltol.*1e-3);
    [t,x] = ode45((@(t,x) EOM(t,x,inputs)),[t1 t2],x1mod,options);
    % [t,x] = ode23((@(t,x) EOM(t,x,inputs)),[t1 t2],x1mod,options);
    % [t,x] = ode113((@(t,x) EOM(t,x,inputs)),[t1 t2],x1mod,options);
%     [t,x] = ode15s((@(t,x) EOM(t,x,inputs)),[t1 t2],x1mod,options);
%     [t,x] = ode23tb((@(t,x) EOM(t,x,inputs)),[t1 t2],x1mod,options);
%     [t,x] = RungeKutta(EOM,[t1 t2],x1mod,inputs.reltol,inputs);
% else                                % Forward run
%     options = odeset('RelTol',inputs.reltol);
    % [t,x] = ode45((@(t,x) EOM(t,x,inputs)),[t1 t2],x1mod,options);        
%     [t,x] = ode23((@(t,x) EOM(t,x,inputs)),[t1 t2],x1mod,options);
%     [t,x] = ode15s((@(t,x) EOM(t,x,inputs)),[t1 t2],x1mod,options);
%     [t,x] = ode113((@(t,x) EOM(t,x,inputs)),[t1 t2],x1mod,options);
% end

x2 = x(end,1:end-covlength)'; % Vector of the propagated state variables
P2vec = x(end,end-covlength+1:end)'; % Vector of the propagated covariance
% Creates the upper-triangle of the eventually symmetric, covariance matrix
P2 = makemat(P2vec);

return