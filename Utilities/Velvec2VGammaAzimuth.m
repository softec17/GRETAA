function [Vmag,gamma,azimuth] = Velvec2VGammaAzimuth(VelVeclist)
% Takes the 3 x 1 velocity vector in the local horizontal coordinate frame 
% and computes
% the velocity magnitude, flight path angle (+ is up) and azimuth angle

num_x = length(VelVeclist);
Vmag = zeros(num_x,1);
gamma = Vmag;
azimuth = Vmag;

for ii = 1:1:num_x

VelVec = VelVeclist(ii,:);
VelVec = VelVec(:);

Vmag(ii) = norm(VelVec);
azimuth(ii) = atan2(VelVec(2),VelVec(1));
gamma(ii) = -atan2(VelVec(3),sqrt(VelVec(1)^2+VelVec(2)^2));
% Recall that + gamma is upwards velocity from the NED frame (as in neg. w)

end

return