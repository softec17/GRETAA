function newacc = AccTimeRate(x,oldacc)

lat   = x(1);                   % rad, geocentric lat
lon = x(2);
q_J2000_2_DS_0 = x(3);          % Quaternions
q_J2000_2_DS_1 = x(4);
q_J2000_2_DS_2 = x(5);
q_J2000_2_DS_3 = x(6);
q_J2000_2_MCMF_0 = x(7);          
q_J2000_2_MCMF_1 = x(8);
q_J2000_2_MCMF_2 = x(9);
q_J2000_2_MCMF_3 = x(10);

q_J2000_DS = [q_J2000_2_DS_0,q_J2000_2_DS_1,q_J2000_2_DS_2,q_J2000_2_DS_3];
q_J2000_MCMF = [q_J2000_2_MCMF_0,q_J2000_2_MCMF_1,q_J2000_2_MCMF_2,q_J2000_2_MCMF_3];
R_J_DS = quat2dcm(q_J2000_DS); R_DS_2_J = R_J_DS';
R_J_MCMF = quat2dcm(q_J2000_MCMF);
R_Lat = [[-sin(lat), 0, cos(lat)]
[                   0, 1,                    0]
[ -cos(lat), 0,  -sin(lat)]];
R_MCMF_LH = R_Lat*R3(lon);
R_DS_LH = R_MCMF_LH*R_J_MCMF*R_DS_2_J; R_LH_DS = R_DS_LH';
newacc = R_LH_DS*oldacc;


return