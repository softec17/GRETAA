function [CA] = aero_interp2_CA(AOA_grid,Mach_grid,data_grid,alpha,beta,Mach)
% Function to interpolate in the aerodynamic database
% alpha, beta = radians
% data_data dimensions in order = CA, CN, Cm, CL, CD

total_alpha = acos(cos(alpha)*cos(beta));
total_alpha = rad2deg(total_alpha);

% Fix the edges
if Mach > max(max(Mach_grid))
    Mach = max(max(Mach_grid));
elseif Mach < min(min(Mach_grid))
    Mach = min(min(Mach_grid));
end

if total_alpha > max(max(AOA_grid))
    total_alpha = max(max(AOA_grid));
elseif total_alpha < min(min(AOA_grid))
    total_alpha = min(min(AOA_grid));
end

CA = interp2(AOA_grid,Mach_grid,data_grid(:,:,1),total_alpha,Mach);

return