function [t,x] = RK4(func,tspan,x0,reltol,varargin)

t = tspan(1):reltol:tspan(end);
x = zeros(length(t),length(x0));
x0 = x0(:);
x(1,:) = x0';

addl_input_flag = 1;
if isempty(varargin)==0
    addl_input_flag = 0;
end

for ii = 2:length(t)
    x_k = x(ii-1,:)';
    t_k = t(ii-1);
    if addl_input_flag == 0
        k1 = func(t_k,x_k);
        k2 = func(t_k + 0.5*reltol,x_k + 0.5*reltol.*k1);
        k3 = func(t_k + 0.5*reltol,x_k + 0.5*reltol.*k2);
        k4 = func(t_k + reltol,x_k + reltol.*k3);
    else
        k1 = func(t_k,x_k,varargin);
        k2 = func(t_k + 0.5*reltol,x_k + 0.5*reltol.*k1,varargin);
        k3 = func(t_k + 0.5*reltol,x_k + 0.5*reltol.*k2,varargin);
        k4 = func(t_k + reltol,x_k + reltol.*k3,varargin);    
    end
    y = (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;
    y = y(:);
    x(ii,:) = y';
end



xdot



return