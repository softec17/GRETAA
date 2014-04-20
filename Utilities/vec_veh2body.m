function [vecrotx,vecroty,vecrotz] = vec_veh2body(vecxin,vecyin,veczin,bankin,...
    pitchin,azimuthin)
% Converts the components of a vector in vehicle frame to body frame
% components
%
% Inputs
% vecxin = x-component of the vector in vehicle frame
% vecyin = y-component of the vector in vehicle frame
% veczin = z-component of the vector in vehicle frame
% azimuthin = azimuth angle (or yaw angle), deg
% pitchin = pitch angle, deg
% bankin = bank angle (or roll angle), deg
%
% Outputs
% vecrotx = x-component of the vector in body frame
% vecroty = y-component of the vector in body frame
% vecrotz = z-component of the vector in body frame

% Length of the vector
N = length(vecxin);

% Preallocate space
vecrotx = zeros(N,1);
vecroty = zeros(N,1);
vecrotz = zeros(N,1);

for ii = 1:1:N
    % Create local copies for the current iteration
    vecx = vecxin(ii);
    vecy = vecyin(ii);
    vecz = veczin(ii);
    bank = bankin(ii);
    pitch = pitchin(ii);
    azimuth = azimuthin(ii);
    % Create a vector out of the components
    vec = [vecx;vecy;vecz];
    % Rotation matrix from vehicle frame to body frame
    Rv2b = R_veh2body(bank,pitch,azimuth);
    % Rotated vector
    vecmod = Rv2b*vec;
    % Split vector into components
    vecrotx(ii) = vecmod(1);
    vecroty(ii) = vecmod(2);
    vecrotz(ii) = vecmod(3);
end

return