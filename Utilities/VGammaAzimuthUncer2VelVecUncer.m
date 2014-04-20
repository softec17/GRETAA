function veluncert = VGammaAzimuthUncer2VelVecUncer(init_V,init_FPA,init_Azimuth,vmaguncer,fpauncer,headinguncer)
% Calculate the uncertainty of the velocity vector in the NED frame if the
% Vmag, gamma, azimuth and its uncertainties are known

velocityVec = VGammaAzimuth2VelVec(init_V,init_FPA,init_Azimuth);
Vn = velocityVec(1);
Ve = velocityVec(2);
Vd = velocityVec(3);

Amat = [Vn, Ve, Vd; Vd*Vn/(Vn^2+Ve^2)^1.5, Vd*Ve/(Vn^2+Ve^2)^1.5, -1/(Vn^2+Ve^2)^0.5; -Ve/Vn^2, 1/Vn, 0];
Bvec = [init_V*vmaguncer;(1/cos(init_FPA)^2)*fpauncer;(1/cos(init_Azimuth)^2)*headinguncer];
veluncert = Amat\Bvec;

return