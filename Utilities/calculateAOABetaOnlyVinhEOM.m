function [AOA,Sideslip] = calculateAOABetaOnlyVinhEOM(state,...
                                         vi,ve,qi,qe,wind_vel)

[num_t,num_c]=size(state);
Vel_veh = state(:,vi:ve);       % Column 1 = Vel, Column 2 = Gamma, Column 3 = Psi(Heading angle)
Quat = state(:,qi:qe);                                     
                                     
% Preallocate space
AOA = zeros(num_t,1);
Sideslip = zeros(num_t,1);
                                     
for ii = 1:1:num_t
    
   % States (AOA and beta)
   q0 = Quat(ii,1); q1 = Quat(ii,2); q2 = Quat(ii,3); q3 = Quat(ii,4);
   R_veh2body = R_veh2body_wq(q0,q1,q2,q3);
   Vel = Vel_veh(ii,1);
   Gamma = Vel_veh(ii,2);
   Psi = Vel_veh(ii,3);
   Vel_vehNED = [Vel*cos(Gamma)*sin(Psi);Vel*cos(Gamma)*cos(Psi);-Vel*sin(Gamma)];
   Vel_body = R_veh2body*Vel_vehNED;
   Vel_wrt_wind = Vel_body - wind_vel(ii,:)';
   [alpha,beta] = Vw_2_alphabeta(Vel_wrt_wind);
   AOA(ii) = deg2rad(alpha);
   Sideslip(ii) = deg2rad(beta);
   
end

return