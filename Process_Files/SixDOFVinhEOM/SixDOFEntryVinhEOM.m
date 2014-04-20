function xdot = SixDOFEntryVinhEOM(t,x,inputs)
% Function calculates the derivative of the states and covaraince at time t
% t = current time (units.time)
% x = current state and covariance matrix elements
%     x(1:3) = position with respect to planet-fixed frame
%     x(4:6) = velocity in vehicle-centered planet-fixed frame (NED)
%     x(7:10) = quaternion componenets (q1,q2,q3,q0)
%     x(11:55) = upper triangle covariance matrix values
%                elements are row by row

n_x = length(x);
n_states = -1.5 + sqrt(2.25+2*n_x);
statedot = SixDOFEntryEOMVinhEOM(t,x(1:n_states),inputs);
covdot = SixDOFEntryCOVVinhEOM(t,x(1+n_states:end),x(1:n_states),inputs);

xdot = [statedot;covdot];

return