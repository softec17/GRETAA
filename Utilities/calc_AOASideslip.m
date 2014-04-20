function [AOA,Sideslip] = calc_AOASideslip(V,gamma,psi,q0,q1,q2,q3)

R_veh2body = R_veh2body_wq(q0,q1,q2,q3);
Vel_vehNED = [V*cos(gamma)*sin(psi);V*cos(gamma)*cos(psi);-V*sin(gamma)];
Vel_body = R_veh2body*Vel_vehNED;
[AOA,Sideslip] = Vw_2_alphabeta(Vel_body);
AOA = deg2rad(AOA);
Sideslip = deg2rad(Sideslip);

return