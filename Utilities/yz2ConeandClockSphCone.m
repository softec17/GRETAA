function [clock,cone] = yz2ConeandClockSphCone(yloc,zloc,rnose,rshell,deltamax)
% Calculate clock and cone angle for pressure port locations for a sphere
% cone
%
% Note: Port locations are assumed to be given based on a body-fixed x =
% out, Z = down right-handed coordinate system
%
% Inputs
% yloc = y coordinates of port (origin at tip of heatshield, axial-out-and
%        down coordinate frame) (meters)
% zloc = z coordinates of port (origin at tip of heatshield, axial-out-and
%        down coordinate frame)  (meters)
% rnose = nose radius (meters)
% rshell = max shell diameter for sphere cone (meters)
% deltamax = cone half angle (deg) 

% % Number of ports
% N = length(yloc);
% if N ~= length(zloc)
%     error('Number of y coordinates of port different than num of z coords.');
% end
% 
% % Limits for cone and sphere
% % Radius limit btw sphere and cone 
% rlim = rnose*cos(deg2rad(deltamax));
% % x distance limit btw sphere and cone from the axis
% xlim = rnose - rnose*sin(deg2rad(deltamax));
% 
% % Radius for every port
% r = sqrt(yloc.^2+zloc.^2);
% 
% % Clock and cone angle
% clock = atan2(-yloc,-zloc);
% cone = zeros(N,1);
% for ii = 1:1:N
%     if r(ii) <= rlim   % If port is on the spherical part
%         cone(ii) = asin(r(ii)/rnose);
%     else               % If port is on the conical part
%         x = xlim + (r(ii) - rlim)/tan(deg2rad(deltamax));
%         cone(ii) = atan((xlim-x*tan(deg2rad(deltamax))-rlim)/(x-rnose));
%     end
% end

% Number of ports
N = length(yloc);
if N ~= length(zloc)
    error('Number of y coordinates of port different than num of z coords.');
end

% Limits for cone and sphere
% Radius limit btw sphere and cone 
rlim = rnose*cos(deg2rad(deltamax));
% x distance limit btw sphere and cone from the axis
xlim = rnose - rnose*sin(deg2rad(deltamax));

% Radius for every port
r = sqrt(yloc.^2+zloc.^2);

% Clock and cone angle
clock = atan2(yloc,zloc);
cone = zeros(N,1);
for ii = 1:1:N
    if r(ii) <= rlim   % If port is on the spherical part
        cone(ii) = asin(r(ii)/rnose);
    else               % If port is on the conical part
        x = xlim + (r(ii) - rlim)/tan(deg2rad(deltamax));
        cone(ii) = atan((xlim-x*tan(deg2rad(deltamax))-rlim)/(x-rnose));
    end
end
cone = abs(asin(sin(cone)));

return

return