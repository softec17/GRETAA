function Uncertainty = processCovariance(Pvec)
% Calculate uncertainties (i.e. 1 sigma values) from covariance matrix
% Assumed that num_col = num_UppertriCovarianceMatrix 
% num_row = num_time_vector
%
% Input
% Pvec = vector of the upper triangle of a covariance matrix
%
% Output
% Uncertainty = vector of 1-sigma uncertainties


[num_t,num_c]=size(Pvec);           % 
num_x = length(makemat(Pvec(1,:))); % Number of states

% Preallocate space
Uncertainty = zeros(num_t,num_x);

for ii = 1:1:num_t
   holder = makemat(Pvec(ii,:));
   holder1 = diag(holder)';
   [i] = find(holder1<0);
   for jj = 1:1:length(i)
       if abs(holder1(i))<1e-10
           holder1(i) = 0;
       end
   end
   holder2 = sqrt(holder1);
   Uncertainty(ii,:) = holder2;
   realcheck = isreal(holder2);
end


return