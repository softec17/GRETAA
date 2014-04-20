function [uq0,uq1,uq2,uq3] = Eul2QuatUncer(bank,pitch,azimuth,Ubank,Upitch,Uazimuth)
% Takes uncertainty in Euler angles and calculates the equivalent 
% uncertainties in the quaternions using chain rule
% Inputs
% bank,pitch,azimuth = Euler angles (rad)in 1, 2, and 3 axis
% ubank,upitch,uazimuth = Euler angles uncertainties (rad)in 1,2,and 3 axis
%
% Outputs
% uq0,uq1,uq2,uq3 = quaternion uncertainties (q0 being the "scalar" value)

theta1 = azimuth;
theta2 = pitch;
theta3 = bank;
utheta1 = Uazimuth;
utheta2 = Upitch;
utheta3 = Ubank;

t1 = theta1/2;
t2 = theta2/2;
t3 = theta3/2;

uq0 = 0.5*(-sin(t1)*cos(t2)*cos(t3)+cos(t1)*sin(t2)*sin(t3))*utheta1+...
      0.5*(-cos(t1)*sin(t2)*cos(t3)+sin(t1)*cos(t2)*sin(t3))*utheta2+...
      0.5*(-cos(t1)*cos(t2)*sin(t3)+sin(t1)*sin(t2)*cos(t3))*utheta3;
uq1 = 0.5*(-sin(t1)*cos(t2)*sin(t3)-cos(t1)*sin(t2)*cos(t3))*utheta1+...
      0.5*(-cos(t1)*sin(t2)*sin(t3)-sin(t1)*cos(t2)*cos(t3))*utheta2+...
      0.5*(cos(t1)*cos(t2)*cos(t3)+sin(t1)*sin(t2)*sin(t3))*utheta3;
uq2 = 0.5*(-sin(t1)*sin(t2)*cos(t3)+cos(t1)*cos(t2)*sin(t3))*utheta1+...
      0.5*(cos(t1)*cos(t2)*cos(t3)-sin(t1)*sin(t2)*sin(t3))*utheta2+...
      0.5*(-cos(t1)*sin(t2)*sin(t3)+sin(t1)*cos(t2)*cos(t3))*utheta3;  
uq3 = 0.5*(cos(t1)*cos(t2)*cos(t3)+sin(t1)*sin(t2)*sin(t3))*utheta1+...
      0.5*(-sin(t1)*sin(t2)*cos(t3)-cos(t1)*cos(t2)*sin(t3))*utheta2+...
      0.5*(-sin(t1)*cos(t2)*sin(t3)-cos(t1)*sin(t2)*cos(t3))*utheta3;    

return