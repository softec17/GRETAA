function [CA,CN,Cm,CL,CD] = aero_interp(aero_data,AOA_list,alpha,beta,Mach)
% Function to interpolate in the aerodynamic database
% alpha, beta = radians
% aero_data columns = Mach, CA, CN, Cm, CL, CD
% 
total_alpha = acos(cos(alpha)*cos(beta));
total_alpha = rad2deg(total_alpha);

[AOA_grid,Mach_grid,new_grid] = create_interp_grid(AOA_list,aero_data);

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

CA = interp2(AOA_grid,Mach_grid,new_grid(:,:,1),total_alpha,Mach,'cubic');
CN = interp2(AOA_grid,Mach_grid,new_grid(:,:,2),total_alpha,Mach,'cubic');
Cm = interp2(AOA_grid,Mach_grid,new_grid(:,:,3),total_alpha,Mach,'cubic');
CL = interp2(AOA_grid,Mach_grid,new_grid(:,:,4),total_alpha,Mach,'cubic');
CD = interp2(AOA_grid,Mach_grid,new_grid(:,:,5),total_alpha,Mach,'cubic');

% AOA_length = length(AOA_list);
% 
% %% Angle of attack look-up
% AOA_ind_lower = find(AOA_list<total_alpha,1,'last');
% AOA_ind_upper = AOA_ind_lower + 1;
% 
% % Check if the AOA is beyond the range of data available
% AOA_lower_empty = 0;
% AOA_upper_empty = 0;
% if isempty(AOA_ind_lower)
%     AOA_lower_empty = 1;
% elseif AOA_ind_lower == AOA_length
%     AOA_upper_empty = 1;
% end
% 
% %% Mach number look-up
% Mach_list = aero_data{1,1}(:,1);        % The Mach number list is the same for all AOAs
% Mach_list_length = length(Mach_list);
% 
% Mach_ind_lower = find(Mach_list<Mach,1,'last');
% Mach_ind_upper = Mach_ind_lower + 1;
% 
% % Check if the mach number is beyond the range of data available
% Mach_lower_empty = 0;
% Mach_upper_empty = 0;
% if isempty(Mach_ind_lower)
%     Mach_lower_empty = 1;
% elseif Mach_ind_lower == Mach_list_length
%     Mach_upper_empty = 1;
% end
% 
% %% Look up aerodynamic values
% if AOA_lower_empty || AOA_upper_empty
%     if AOA_lower_empty                  % AOA is lower than available data -- use lowest available data
%         if Mach_lower_empty || Mach_upper_empty
%             if Mach_lower_empty         % Mach number is lower than available data -- use lowest available data
%                 list_aero = aero_data{1,1}(1,2:end);
%             else                        % Mach number is higher than available data -- use highest available data
%                 list_aero = aero_data{1,1}(Mach_list_length,2:end);
%             end
%         else
%             Mach_scale_factor = (Mach - Mach_list(Mach_ind_lower))/(Mach_list(Mach_ind_upper) - Mach_list(Mach_ind_lower));
%             list_lower_Mach = aero_data{1,1}(Mach_ind_lower,2:end);
%             list_upper_Mach = aero_data{1,1}(Mach_ind_upper,2:end);
%             list_aero = Mach_scale_factor.*(list_upper_Mach - list_lower_Mach) + list_lower_Mach;
%         end
%     else                                % AOA is higher than available data -- use highest available data
%         if Mach_lower_empty || Mach_upper_empty
%             if Mach_lower_empty         % Mach number is lower than available data -- use lowest available data
%                 list_aero = aero_data{AOA_length,1}(1,2:end);
%             else                        % Mach number is higher than available data -- use highest available data
%                 list_aero = aero_data{AOA_length,1}(Mach_list_length,2:end);
%             end
%         else
%             Mach_scale_factor = (Mach - Mach_list(Mach_ind_lower))/(Mach_list(Mach_ind_upper) - Mach_list(Mach_ind_lower));
%             list_lower_Mach = aero_data{AOA_length,1}(Mach_ind_lower,2:end);
%             list_upper_Mach = aero_data{AOA_length,1}(Mach_ind_upper,2:end);
%             list_aero = Mach_scale_factor.*(list_upper_Mach - list_lower_Mach) + list_lower_Mach;
%         end
%     end
% else
%     if Mach_lower_empty || Mach_upper_empty
%         if Mach_lower_empty             % Mach number is lower than available data -- use lowest available data
%             AOA_scale_factor = (total_alpha - AOA_list(AOA_ind_lower))/(AOA_list(AOA_ind_upper) - AOA_list(AOA_ind_lower));
%             list_lower_AOA = aero_data{AOA_ind_lower,1}(1,2:end);
%             list_upper_AOA = aero_data{AOA_ind_upper,1}(1,2:end);
%             list_aero = AOA_scale_factor.*(list_upper_AOA - list_lower_AOA) + list_lower_AOA;
%         else                            % Mach number is higher than available data -- use highest available data
%             AOA_scale_factor = (total_alpha - AOA_list(AOA_ind_lower))/(AOA_list(AOA_ind_upper) - AOA_list(AOA_ind_lower));
%             list_lower_AOA = aero_data{AOA_ind_lower,1}(Mach_list_length,2:end);
%             list_upper_AOA = aero_data{AOA_ind_upper,1}(Mach_list_length,2:end);
%             list_aero = AOA_scale_factor.*(list_upper_AOA - list_lower_AOA) + list_lower_AOA;
%         end
%     else
%         AOA_scale_factor = (total_alpha - AOA_list(AOA_ind_lower))/(AOA_list(AOA_ind_upper) - AOA_list(AOA_ind_lower));
%         Mach_scale_factor = (Mach - Mach_list(Mach_ind_lower))/(Mach_list(Mach_ind_upper) - Mach_list(Mach_ind_lower));
%         
%         list_lower_AOA_lower_Mach = aero_data{AOA_ind_lower,1}(Mach_ind_lower,2:end);
%         list_lower_AOA_upper_Mach = aero_data{AOA_ind_lower,1}(Mach_ind_upper,2:end);
%         list_lower_AOA = Mach_scale_factor.*(list_lower_AOA_upper_Mach - list_lower_AOA_lower_Mach) + list_lower_AOA_lower_Mach;
%         
%         list_upper_AOA_lower_Mach = aero_data{AOA_ind_upper,1}(Mach_ind_lower,2:end);
%         list_upper_AOA_upper_Mach = aero_data{AOA_ind_upper,1}(Mach_ind_upper,2:end);
%         list_upper_AOA = Mach_scale_factor.*(list_upper_AOA_upper_Mach - list_upper_AOA_lower_Mach) + list_upper_AOA_lower_Mach;
%         
%         list_aero = AOA_scale_factor.*(list_upper_AOA - list_lower_AOA) + list_lower_AOA;
%     end
% end
% 
% CA = list_aero(1);
% CN = list_aero(2);
% Cm = list_aero(3);
% CL = list_aero(4);
% CD = list_aero(5);

return