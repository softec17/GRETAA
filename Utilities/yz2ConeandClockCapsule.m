function [clock,cone] = yz2ConeandClockCapsule(yloc,zloc,rnose,deltamax)
% Calculate clock and cone angle for pressure port locations
%
% Inputs
% yloc = y coordinates of port (origin at tip of heatshield, axial-out-and
%        down coordinate frame)
% zloc = z coordinates of port (origin at tip of heatshield, axial-out-and
%        down coordinate frame)
% rnose = nose radius
% deltamax = cone half angle (deg.)

% Number of ports
N = length(yloc);
if N ~= length(zloc)
    error('Number of y coordinates of port different than num of z coords.');
end

% Radius
r = norm_vec([yloc,zloc]);

% Clock and cone angle
clock = atan2(zloc,yloc);
cone = asin(r./rnose);

% Check for errors
ind = find(cone>deg2rad(deltamax),1);
if isempty(ind) ~= 1
    keyboard
    error('Port locations not on front shell of the heatshield');
end

return