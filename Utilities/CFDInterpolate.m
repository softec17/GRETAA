function [Cp] = CFDInterpolate(Clock,Cone,alpha,beta,Mach,aero_data)
% Looks up the coefficient of pressure at the pressure port location
% Inputs
% Clock = clock angle (rad)
% Cone = cone angle (rad)
% alpha = angle-of-attack (rad)
% beta = side-slip angle (rad)
% Mach = mach number of the freestream
% aero_data = Cp tables from CFD
% Output
% Cp = coefficient of pressure
%
% This file is heavily adapted from AMA's newStep code

% Wind relative velocity direction unit vector
Wind_x = cos(alpha)*cos(beta);
Wind_y = sin(beta);
Wind_z = sin(alpha)*cos(beta);

% Total angle of attack
tot_AOA = atan2(sqrt(Wind_y^2+Wind_z^2),Wind_x);

% Wind axis bank angle
bank = atan2(Wind_y,Wind_z);

% Subtract the effect of the bank angle from clock angle
Clock = Clock - bank;

% Restrict Clock angle between -pi and pi
% Adapted from AMA's newStep code
Clock(Clock>pi) = Clock(Clock>pi) - 2*pi;
Clock(Clock<-pi) = Clock(Clock<-pi) + 2*pi;
% Reflect negative clock angles positive angles
% Body is axysymeetric and alpha effects can be explained by tot_AOA in the
% look-up tables
Clock = abs(Clock);

% Mach number vector from the aerodynamics dataset
Machlist = aero_data.mach_vec;
% Clock angle vector from the aerodynamics dataset
Clocklist = aero_data.clock_vec';
% Cone angle vector from the aerodynamics dataset
Conelist = aero_data.cone_vec';

% Convert all angles to degrees (table is in degrees)
tot_AOA = rad2deg(tot_AOA);
Clock = rad2deg(Clock);
Cone = rad2deg(Cone);

% Bounding indices
M_1ind = find(Machlist<=Mach,1,'last');
M_1ind(Mach>=Machlist(end)) = size(Machlist,1)-1;
M_1ind(Mach<Machlist(1)) = 1;
M_2ind = M_1ind + 1;
% M_2ind = find(Machlist>=Mach,1,'first');
% if isempty(M_2ind)
%     M_2ind = M_1ind;
% end
num_ports = length(Clock);
Cp = zeros(1,num_ports);

for jj = 1:1:num_ports
    if Machlist(M_1ind)~=Mach
        AOAlist = aero_data.mach(M_2ind).aoa;
        Cplist = aero_data.mach(M_2ind).Cp;
        dim = size(Cplist);
        if ((dim(1)~=length(AOAlist))||(dim(2)~=length(Clocklist))||(dim(3)~=length(Conelist)))
            keyboard
        else
            Cp(jj) = sinterp3(AOAlist,Clocklist,Conelist,Cplist,tot_AOA,Clock(jj),...
                Cone(jj));
        end
    else
        AOAlist = aero_data.mach(M_1ind).aoa;
        Cplist = aero_data.mach(M_1ind).Cp;
        Cp(jj) = sinterp3(AOAlist,Clocklist,Conelist,Cplist,tot_AOA,Clock(jj),...
            Cone(jj));
    end
end

return