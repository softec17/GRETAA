function [yaw,pitch,roll] = InerQuat2Eul(q0,q1,q2,q3,lat,lon,PrimeMeridian)
% Takes quaternions (inertial to body) and calculates the equivalent Euler angles
% Inputs
% q0,q1,q2,q3 = quaternions with q0 being the "scalar" value
%
% Outputs
% yaw,pitch,roll = Euler angles (rad) in 3,2,and 1 axis

q0 = q0(:);
q1 = q1(:);
q2 = q2(:);
q3 = q3(:);
lat = lat(:);
lon = lon(:);
PrimeMeridian = PrimeMeridian(:);

yaw = zeros(size(q0));
pitch = zeros(size(q0));
roll = zeros(size(q0));

for ii = 1:1:length(yaw)
    RIner2Cru = quat2dcm([q0(ii),q1(ii),q2(ii),q3(ii)]);
    RIner2Vehmat = R_iner2veh(lat(ii),lon(ii)+PrimeMeridian(ii));
    RVeh2Cru = RIner2Vehmat'*RIner2Cru;
    RCru2Bdy = [1 0 0;0 0 -1;0 1 0];
    RVeh2Bdy = RCru2Bdy*RVeh2Cru;
    [yawh, pitchh, rollh] = dcm2angle(RVeh2Bdy);
    yaw(ii) = yawh;
    pitch(ii) = pitchh;
    roll(ii) = rollh;
end

return

