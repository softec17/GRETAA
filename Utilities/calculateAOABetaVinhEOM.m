function [AOA,Sideslip,UAOA,USideslip] = calculateAOABetaVinhEOM(state,P_vec,...
                                         vi,ve,qi,qe,wind_vel,lat,lon,primemeridian)
% Calculate AOA and beta and their uncertainties
% Assumed that num_col = num_UppertriCovarianceMatrix 
% num_row = num_time_vector
% Assumed wind velocity is expressed in the body frame
%
% Inputs
% state = state vector [position;velocity;quaternions;pressure;density]
% P_vec = vector representation of upper-tri of covariance matrix
% vi, ve = index for start and end of velocity in state vector
% qi, qe = index for start and end of quaternions in state vector
% wind_vel = wind speed in body frame
% lat = latitude
% lon = longitude
% primemeridian = prime meriadian's angle from the right ascension
% Outputs
% AOA = angle of attack
% Sideslip = sideslip angle
% UAOA = uncertainty in the angle of attack
% USideslip = uncertainty in sideslip angle

[num_t,num_c]=size(state);

% Preallocate space
AOA = zeros(num_t,1);
Sideslip = zeros(num_t,1);
UAOA = zeros(num_t,1);
USideslip = zeros(num_t,1);

Vel_veh = state(:,vi:ve);       % Column 1 = Vel, Column 2 = Gamma, Column 3 = Psi(Heading angle)
Quat = state(:,qi:qe);

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
   
   % Uncertainty (in AOA and beta)
   holder = makemat(P_vec(ii,:));
   Pvel = holder(vi:ve,vi:ve);
   Pquat = holder(qi:qe,qi:qe);
   
   holder2 = diag(Pvel);
   realcheck = isreal(holder2);
   if realcheck ~= 1
       keyboard 
       error(['Imaginary standard deviation found. ii = ',num2str(ii)])
   end
   Uvel = holder(4,4);
   Uvelgamma = holder(5,5);
   Uvelpsi = holder(6,6);   %sqrt(holder2);
   
   % NED velocity uncertainties
   % NED velocities
    u = Vel*cos(Gamma)*sin(Psi);
    v = Vel*cos(Gamma)*cos(Psi);
    w = -Vel*sin(Gamma);
   Uvelveh = zeros(3,1);
   Uvelveh(1) = cos(Gamma)*sin(Psi)*Uvel + Vel*cos(Gamma)*cos(Psi)*Uvelpsi - Vel*sin(Gamma)*sin(Psi)*Uvelgamma;
   Uvelveh(2) = cos(Gamma)*cos(Psi)*Uvel - Vel*cos(Gamma)*sin(Psi)*Uvelpsi - Vel*sin(Gamma)*cos(Psi)*Uvelgamma;
   Uvelveh(3) = -sin(Gamma)*Uvel - Vel*cos(Gamma)*Uvelgamma;
   
   holder3 = diag(Pquat);
   realcheck = isreal(holder3);
   if realcheck ~= 1
       keyboard 
       error(['Imaginary standard deviation found. ii = ',num2str(ii)])
   end
   Uquat = sqrt(holder3);
   
   dq0 = [4*q0,2*q3,-2*q2;-2*q3,4*q0,2*q1;2*q2,-2*q1,4*q0]*Vel_veh(ii,:)';
   dq1 = [4*q1,2*q2,2*q3;2*q2,0,2*q0;2*q1,-2*q0,0]*Vel_veh(ii,:)';
   dq2 = [0,2*q1,-2*q0;2*q1,4*q2,2*q3;2*q0,2*q3,0]*Vel_veh(ii,:)';
   dq3 = [0,2*q0,2*q1;-2*q0,0,2*q2;2*q1,2*q2,4*q3]*Vel_veh(ii,:)';
   dVvehx = [2*q0^2-1+2*q1^2;2*(q1*q2-q0*q3);2*(q1*q3+q0*q2)];
   dVvehy = [2*(q1*q2+q0*q3);2*q0^2-1+2*q2^2;2*(q2*q3-q0*q1)];
   dVvehz = [2*(q1*q3-q0*q2);2*(q2*q3+q0*q1);2*q0^2-1+2*q3^2];

   Uvelb = abs(dq0*Uquat(1)+dq1*Uquat(2)+dq2*Uquat(3)+dq3*Uquat(4)+dVvehx*Uvelveh(1)+dVvehy*Uvelveh(2)+dVvehz*Uvelveh(3));
%    Uvelb = abs(dq0*Uquat(1)+dq1*Uquat(2)+dq2*Uquat(3)+dq3*Uquat(4)+Uvelveh);
   Uvelbnorm = norm(Uvelb);
   
   V_inf = Vel_wrt_wind;
   UAOA(ii) = (1/sec(AOA(ii))^2)*abs(Uvelb(3)/V_inf(1)+(V_inf(3)/V_inf(1)^2)*Uvelb(1));
   USideslip(ii) = abs((1/cos(Sideslip(ii)))*(Uvelb(2)/norm(V_inf)+(V_inf(2)/norm(V_inf)^2)*Uvelbnorm));
% 
%     sigu = [(2 * q0 ^ 2 + 2 * q1 ^ 2 - 1) * cos(gmma) * sin(psi) + (2 * q1 * q2 + 2 * q0 * q3) * cos(gmma) * cos(psi) - (2 * q1 * q3 - 2 * q0 * q2) * sin(gmma) -(2 * q0 ^ 2 + 2 * q1 ^ 2 - 1) * V * sin(gmma) * sin(psi) - (2 * q1 * q2 + 2 * q0 * q3) * V * sin(gmma) * cos(psi) - (2 * q1 * q3 - 2 * q0 * q2) * V * cos(gmma) (2 * q0 ^ 2 + 2 * q1 ^ 2 - 1) * V * cos(gmma) * cos(psi) - (2 * q1 * q2 + 2 * q0 * q3) * V * cos(gmma) * sin(psi) 0.4e1 * q0 * V * cos(gmma) * sin(psi) + 0.2e1 * q3 * V * cos(gmma) * cos(psi) + 0.2e1 * q2 * V * sin(gmma) 0.4e1 * q1 * V * cos(gmma) * sin(psi) + 0.2e1 * q2 * V * cos(gmma) * cos(psi) - 0.2e1 * q3 * V * sin(gmma) 0.2e1 * q1 * V * cos(gmma) * cos(psi) + 0.2e1 * q0 * V * sin(gmma) 0.2e1 * q0 * V * cos(gmma) * cos(psi) - 0.2e1 * q1 * V * sin(gmma)]*...
%         [sigV;sigGamma;sigPsi;sigq0;sig21;sigq2;sigq3];




      
end

return