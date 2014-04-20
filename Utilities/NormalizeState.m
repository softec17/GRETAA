function [x] = NormalizeState(xold,inputs)

x = xold;

% normalize.r = 1.27878107345943E+05 + radius_planet;
% normalize.V = 6000;
% normalize.dens = 0.02;
% normalize.press = 636;

normalize = inputs.normalize;

x(1) = x(1)./normalize.r;
x(4) = x(4)./normalize.V;
x(11) = x(11)./normalize.press;
x(12) = x(12)./normalize.dens;

return