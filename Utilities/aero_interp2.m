function [CA,CN,Cm,CL,CD] = aero_interp2(AOA_grid,Mach_grid,data_grid,alpha,beta,Mach)
% Function to interpolate in the aerodynamic database
% alpha, beta = radians
% data_data dimensions in order = CA, CN, Cm, CL, CD

total_alpha = acos(cos(alpha)*cos(beta));
total_alpha = rad2deg(total_alpha);

% Fix the edges
if Mach > max(max(Mach_grid))
    Mach = max(max(Mach_grid));
%     keyboard
elseif Mach < min(min(Mach_grid))
    Mach = min(min(Mach_grid));
%     keyboard
end

if total_alpha > max(max(AOA_grid))
    total_alpha = max(max(AOA_grid));
%     keyboard
elseif total_alpha < min(min(AOA_grid))
    total_alpha = min(min(AOA_grid));
%     keyboard
end

CA = interp2(AOA_grid,Mach_grid,data_grid(:,:,1),total_alpha,Mach);
CN = interp2(AOA_grid,Mach_grid,data_grid(:,:,2),total_alpha,Mach);
Cm = interp2(AOA_grid,Mach_grid,data_grid(:,:,3),total_alpha,Mach);
CL = interp2(AOA_grid,Mach_grid,data_grid(:,:,4),total_alpha,Mach);
CD = interp2(AOA_grid,Mach_grid,data_grid(:,:,5),total_alpha,Mach);

return