function xsig = CreateSigmaPt(x,P,kappa)
% Creates (2n-1) sigma points based on nx1 state vector x and nxn
% covariance matrix. kappa is an user-defined tuning parameter
% See notation in Julier et al. 2000 IEEE Transactions document 

n = length(x);
xsig = zeros(n,2*n+1);
alpha = 1e-3;
lamda = alpha^2*(n+kappa)-n;
pert_P = (n+lamda).*P;
sqrt_pert_P_trans = chol(pert_P);
% sqrt_pert_P_trans = sqrtm(pert_P);
if isreal(sqrt_pert_P_trans)~=1
    keyboard;
end
sqrt_pert_P_trans = real(sqrt_pert_P_trans);

xsig(:,1) = x;
for ii = 1:1:n
   xsig(:,1+ii) = x + sqrt_pert_P_trans(ii,:)';
   xsig(:,1+n+ii) = x - sqrt_pert_P_trans(ii,:)';
end
if isreal(xsig)~=1
    keyboard;
end

return