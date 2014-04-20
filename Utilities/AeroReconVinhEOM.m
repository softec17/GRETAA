%% Aerodynamic reconstruction using IMU data

% Rate of reconstruction
raterecon = 1e-3;

% Mass time line
aerorecontime = fwdstate.time(1):raterecon:fwdstate.time(end);
aeroreconmass = interp1(actual.time(indices),actual.mass(indices),aerorecontime)';

% Reference Area (assumed to be base area)
aeroreconarea = (pi/4)*2.65^2;

% Dynamic pressure reconstructed
aerorecondens = interp1(fwdstate.time,fwdstate.state(:,12),aerorecontime)';
aeroreconvel = interp1(fwdstate.time,fwdstate.state(:,4),aerorecontime)';
aeroreconq = 0.5.*aerorecondens.*aeroreconvel.^2;
aerorecondensuncer = interp1(fwdstate.time,Uncervec(:,12),aerorecontime)';
aeroreconveluncer = interp1(fwdstate.time,Uncervec(:,4),aerorecontime)';

% Acceleration time history
aeroreconaccel = interp1(timeAccel,AccelSampled,aerorecontime)';

%% Axial Force
% Axial force is from IMU data in the negative x axis (body frame)
AxialForce = aeroreconmass.*aeroreconaccel(:,1);

% Reconstructed force coefficients
CAEst = AxialForce./(aeroreconq.*aeroreconarea);

% Reconstruct force coefficient uncertainties
CAUncer = ((2.*aeroreconmass)./(aeroreconarea.*aeroreconvel.^2.*aerorecondens)).*IMU_accel_uncertainty + ...
    -((4.*aeroreconaccel.*aeroreconmass)./(aeroreconarea.*aeroreconvel.^3.*aerorecondens)).*aeroreconvel + ...
    -((2.*aeroreconaccel.*aeroreconmass)./(aeroreconarea.*aeroreconvel.^2.*aerorecondens.^2)).*aerorecondensuncer;

% Axial Force Coefficient Reconstruction plot
figure(60)
clf
set(60,'WindowStyle','docked','name','CA');
subplot(2,1,1)
plot(actual.time(indices),actual.ca(indices),'r','LineWidth',2)
hold on
plot(aerorecontime(1:1000:end),CAEst(1:1000:end),'ko','MarkerSize',6);
hold off
xlabel('Time (sec)','FontWeight','bold')
ylabel('Axial force coefficient, C_A (unitless)','FontWeight','bold')
grid on
legend('Actual C_A','Recon C_A');
title('CA vs. time','FontSize',fontsize,'FontWeight',fontweight)
subplot(2,1,2)
actualca = interp1(actual.time(indices),actual.ca(indices),aerorecontime)';
diffca = actualca - CAEst;
plot(aerorecontime,diffca,'k--','LineWidth',2)
hold on
plot(aerorecontime,CAUncer(:,1),'g--','LineWidth',2);
plot(aerorecontime,-CAUncer(:,1),'g--','LineWidth',2);
plot(aerorecontime,diffca,'k--','LineWidth',2)
% axis([0 300 -1 1]);
hold off
legend('Residual','Est. \pm1\sigma uncertainty','Location','Best')
xlabel('Time (sec)','FontSize',fontsize,'FontWeight',fontweight)
ylabel('Residual CA (unitless)','FontSize',fontsize,'FontWeight',fontweight)
title('Residual CA vs. time','FontSize',fontsize,'FontWeight',fontweight)
grid on

%% Normal Force
% Normal force is from IMU data in the y and z axis actual
% CFD frame is different than body frame for trajectory
% CFD frame in terms of body frame: x = -x, y = y, z = -z
aeroreconaccel = aeroreconaccel';
NormalForce = aeroreconmass.*(-1).*sign(aeroreconaccel(:,3)).*sqrt(aeroreconaccel(:,3).^2);

% Reconstructed force coefficients
CNEst = NormalForce./(aeroreconq.*aeroreconarea);

% Reconstruct force coefficient uncertainties
CNUncer = (aeroreconmass./(aeroreconq.*aeroreconarea)).*IMU_accel_uncertainty;% + ...
%             (NormalForce./aeroreconarea).*(-1./aeroreconq.^2).*aeroreconquncer;

%% Normal Force Coefficient Reconstruction
figure(61)
clf
set(61,'WindowStyle','docked','name','CN');
subplot(2,1,1)
plot(actual.time(indices),actual.cn(indices),'r','LineWidth',2)
% plot(actual.time,actual.cn_nom,'b--','LineWidth',2)
hold on
plot(aerorecontime(1:1000:end),CNEst(1:1000:end),'ko','MarkerSize',6);
hold off
xlabel('Time (sec)','FontWeight','bold')
ylabel('Normal force coefficient, C_N (unitless)','FontWeight','bold')
grid on
legend('Actual C_N','Recon C_N');
title('CN vs. time','FontSize',fontsize,'FontWeight',fontweight)
subplot(2,1,2)
actualcn = interp1(actual.time(indices),actual.cn(indices),aerorecontime)';
diffcn = actualcn - CNEst;
plot(aerorecontime,diffcn,'k--','LineWidth',2)
hold on
plot(aerorecontime,CNUncer(:,1),'g--','LineWidth',2);
plot(aerorecontime,-CNUncer(:,1),'g--','LineWidth',2);
plot(aerorecontime,diffcn,'k--','LineWidth',2)
% axis([0 300 -1 1]);
hold off
legend('Residual','Est. \pm1\sigma uncertainty','Location','Best')
xlabel('Time (sec)','FontSize',fontsize,'FontWeight',fontweight)
ylabel('Residual CN (unitless)','FontSize',fontsize,'FontWeight',fontweight)
title('Residual CN vs. time','FontSize',fontsize,'FontWeight',fontweight)
grid on

