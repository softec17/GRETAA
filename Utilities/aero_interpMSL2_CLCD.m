function [CL,CD] = aero_interpMSL2_CLCD(AOA_grid,Mach_grid,data_grid,alpha,beta,Mach)
% Function to interpolate in the aerodynamic database
% alpha, beta = radians
% data_data dimensions in order = CA, CN, Cm

total_alpha = acos(cos(alpha)*cos(beta));
total_alpha = rad2deg(total_alpha);

if Mach >= 10.63                            % Hypersonic
    AOA_grid_local = AOA_grid.AOA_grid_hyp;
    Mach_grid_local = Mach_grid.Mach_grid_hyp;
    data_grid_local = data_grid.data_grid_hyp;
elseif Mach > 8.86 && Mach < 10.63          % Transitional
    data_grid_local = [data_grid.data_grid_8_86;data_grid.data_grid_hyp(1,[1:3,5:6],:)];
    AOA_grid_local = [AOA_grid.AOA_list_8_86';AOA_grid.AOA_list_8_86'];
    Mach_grid_local = [Mach_grid.Mach_grid_hyp(1,[1:3,5:6]);8.86.*ones(1,5)];
elseif Mach > 5.83 && Mach < 8.86           % Transitional
    data_grid_local = [data_grid.data_grid_sub(5,4:6,:);data_grid.data_grid_8_86(1,[2,4:5],:)];
    AOA_grid_local = [AOA_grid.AOA_list_8_86([2,4:5])';AOA_grid.AOA_list_8_86([2,4:5])'];
    Mach_grid_local = [Mach_grid.Mach_grid_hyp(5,4:6);8.86.*ones(1,3)];
else                                        % Supersonic and subsonic
    AOA_grid_local = AOA_grid.AOA_grid_sub;
    Mach_grid_local = Mach_grid.Mach_grid_sub;
    data_grid_local = data_grid.data_grid_sub;
end

% Fix the edges for Mach number
if Mach > max(max(Mach_grid_local))
    Mach = max(max(Mach_grid_local));
elseif Mach < min(min(Mach_grid_local))
    Mach = min(min(Mach_grid_local));
end

% Fix the edges for angle of attack
if total_alpha > max(max(AOA_grid_local))
    total_alpha = max(max(AOA_grid_local));
elseif total_alpha < min(min(AOA_grid_local))
    total_alpha = min(min(AOA_grid_local));
end

[CL,CD] = aero_interp2_small(AOA_grid_local,Mach_grid_local,data_grid_local,total_alpha,Mach);

return

function [CL,CD] = aero_interp2_small(AOA_grid,Mach_grid,data_grid,total_alpha,Mach)
CA = interp2(AOA_grid,Mach_grid,data_grid(:,:,1),total_alpha,Mach);
CN = interp2(AOA_grid,Mach_grid,data_grid(:,:,2),total_alpha,Mach);
CL = -CA*sind(total_alpha) + CN*cosd(total_alpha);
CD = CA*cosd(total_alpha) + CN*sind(total_alpha);
if CD < 0
    keyboard
end
return