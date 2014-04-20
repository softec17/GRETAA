function [t,x] = RungeKutta(EOM,tspan,x0,reltol,inputs)
% 4th Order Runge-Kutta

t1 = tspan(1);
t2 = tspan(2);
t = [t1:reltol:t2]';
lent = length(t);
x = zeros(lent,length(x0)); x(1,:) = x0';

for ii = 2:lent
    x1 = x(ii-1,:)';
    k1 = EOM(t(ii-1),x1,inputs);
    k2 = EOM(t(ii-1)+0.5*reltol,x1 + 0.5.*reltol.*k1,inputs);
    k3 = EOM(t(ii-1)+0.5*reltol,x1 + 0.5.*reltol.*k2,inputs);
    k4 = EOM(t(ii),x1 + reltol.*k3,inputs);
    x2 = x1 + reltol.*((1/6)*(k1 + k4) + (1/3)*(k2 + k3));
    x(ii,:) = x2';
end

return