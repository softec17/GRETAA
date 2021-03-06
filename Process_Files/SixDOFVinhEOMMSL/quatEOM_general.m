function xdot =  quatEOM_general(t,x,omega)

q0 = x(1); q1 = x(2); q2 = x(3); q3 = x(4);

xdot = 0.5*[-q1 -q2 -q3;...
    q0 -q3 q2;...
    q3 q0 -q1;...
    -q2 q1 q0;]*omega;

return