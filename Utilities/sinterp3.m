function [vi] = sinterp3(x,y,z,v,xi,yi,zi)
% This code is the exact 3D look-up code developed by AMA and used by them
% extensively for Cp lookup tables to generate MEADS-like data
% This function was provided by AMA when Cp tables were provided
% FUNCTION sinterp3
% Simple function to do 3-D linear interpolation
% ASSUMES that x,y,z are monotonic increasing vectors, 
%   v is (length(x),length(y),length(z)) in size,
%   xi,yi,zi are vectors of same length, 
%   performs linear interpolation and extrapolation,
%   no error checking performed for speed reasons
% OUTPUTS vi is a column vector of length(xi)
%
% *** WARNING ***
% This code has NO ERROR CHECKING AT ALL, it was built for speed

% Get indices of x,y,z just less than (x1) and greater than (x2) xi,yi,zi
x = x(:);
xi = xi(:);
[ignore,x1] = histc(xi,x);
if any(x1==0)||any(x1==length(x))
    x1(xi>=x(end)) = size(x,1)-1;
    x1(xi<x(1)) = 1;
end
x2 = x1 + 1;
nrows = length(x);

y = y(:);
yi = yi(:);
[ignore,y1] = histc(yi,y);
if any(y1==0)||any(y1==length(y))
    y1(yi>=y(end)) = size(y,1)-1;
    y1(yi<y(1)) = 1;
end
y2 = y1 + 1;
ncols = length(y);

z = z(:);
zi = zi(:);
[ignore,z1] = histc(zi,z);
if any(z1==0)||any(z1==length(z))
    z1(zi>=z(end)) = size(z,1)-1;
    z1(zi<z(1)) = 1;
end
z2 = z1 + 1;

%For single inputs (xi,yi,zi):
%Interpolate on z
% rz = (zi-z(z1))./(z(z2)-z(z1));
% iz11 = v(x1,y1,z1)+rz.*(v(x1,y1,z2)-v(x1,y1,z1));
% iz12 = v(x1,y2,z1)+rz.*(v(x1,y2,z2)-v(x1,y2,z1));
% iz21 = v(x2,y1,z1)+rz.*(v(x2,y1,z2)-v(x2,y1,z1));
% iz22 = v(x2,y2,z1)+rz.*(v(x2,y2,z2)-v(x2,y2,z1));
%Interpolate on y
% ry = (yi-y(y1))./(y(y2)-y(y1));
% iy1 = iz11+ry.*(iz12-iz11);
% iy2 = iz21+ry.*(iz22-iz21);
%Interpolate on x
% rx = (xi-x(x1))./(x(x2)-x(x1));
% ix = iy1+rx.*(iy2-iy1);
% vi = ix;

%For multiple (xi,yi,zi):
%Linear indexing formula derived from sub2ind.m for 3 dimensions
% ((1+(x1-1)*1)+(y1-1)*nrows)+(z1-1)*nrows*ncols and simplified
%x and y equations simplified to one
rz = (zi-z(z1))./(z(z2)-z(z1));
ry = (yi-y(y1))./(y(y2)-y(y1));
iz111ind = x1+y1*nrows+z1*nrows*ncols-nrows*(1+ncols);
iz11 = v(iz111ind)+rz.*(v(x1+y1*nrows+z2*nrows*ncols-nrows*(1+ncols))-v(iz111ind));
iz121ind = x1+y2*nrows+z1*nrows*ncols-nrows*(1+ncols);
iz12 = v(iz121ind)+rz.*(v(x1+y2*nrows+z2*nrows*ncols-nrows*(1+ncols))-v(iz121ind));
iz211ind = x2+y1*nrows+z1*nrows*ncols-nrows*(1+ncols);
iz21 = v(iz211ind)+rz.*(v(x2+y1*nrows+z2*nrows*ncols-nrows*(1+ncols))-v(iz211ind));
iz221ind = x2+y2*nrows+z1*nrows*ncols-nrows*(1+ncols);
iz22 = v(iz221ind)+rz.*(v(x2+y2*nrows+z2*nrows*ncols-nrows*(1+ncols))-v(iz221ind));
vi = (iz11+ry.*(iz12-iz11))+(xi-x(x1))./(x(x2)-x(x1)).*((iz21+ry.*(iz22-iz21))-(iz11+ry.*(iz12-iz11)));
