function velocityVec = VGammaAzimuth2VelVec(V,Gamma,Azimuth)
% Angles are assumed in radians
% Takes a velocity magnitude (V), inertial flight path angle (gamma) and
% heading angle (azimuth) and converts it into a velocity vector

% Note: Flight path angle convention is neg. for downward traj (re-entry)
% and pos. for upward traj (lifting). But the EOM coordinate system is
% North-East-Down, so the sign of the flight path angle must be reversed
Gamma = Gamma*-1;

velocityVec = zeros(3,1);

velocityVec(1) = V*cos(Gamma)*cos(Azimuth);
velocityVec(2) = V*cos(Gamma)*sin(Azimuth);
velocityVec(3) = V*sin(Gamma);

return