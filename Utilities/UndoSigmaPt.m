function [x,P] = UndoSigmaPt(xsig,kappa)
% Creates the mean (x) and the covariance matrix (P) from the sigma pts 
% (xsig) and uses kappa which is an user defined weighting function

% # of states in x and P
% Assume that # of rows of xsig is n
n = size(xsig,1);
alpha = 1e-3;
beta = 2;
lamda = alpha^2*(n+kappa)-n;
x = zeros(n,1);
P = zeros(n);
W_0_mean = lamda/(n+lamda);
W_0_cov = W_0_mean + 1 - alpha^2 + beta;
W_i = 1/(2*(n+lamda));

for ii = 1:1:2*n+1
   if ii == 1
       x = x + W_0_mean.*xsig(:,ii);
   else
       x = x + W_i.*xsig(:,ii);
   end
end
for ii = 1:1:2*n+1
   if ii == 1
      var = (xsig(:,ii)-x); 
      P = P + W_0_cov.*(var*var');
   else
      var = (xsig(:,ii)-x); 
      P = P + W_i.*(var*var');
   end
end

return