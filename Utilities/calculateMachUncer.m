function Machuncer = calculateMachUncer(Pvec,Vvec,i,e)
% Calculate the uncertainty in the Mach number
% Assumed that num_col = num_UppertriCovarianceMatrix 
% num_row = num_time_vector
%
% Inputs
% Pvec = vector representation of covariance matrix 
% Vvec = velocity vector
% i, e = index of start and stop of the velocities in the covariance matrix
% Output
% Machuncer = uncertainty in the Mach number

[num_t,num_c] = size(Pvec);

% Preallocate space
Machuncer = zeros(num_t,1);

for ii = 1:1:num_t
   Mach = norm([Vvec(ii,1),Vvec(ii,2),Vvec(ii,3)]);
   holder = makemat(Pvec(ii,:));
   holder1 = holder(i:e,i:e);
   holder2 = eig(holder1);
   realcheck = isreal(holder2);
   if realcheck ~= 1
       keyboard 
       error(['Imaginary standard deviation found. ii = ',num2str(ii)])
   end
   holder2 = sqrt(holder2);     % Turn covariance into standard deviation
   Machuncer(ii) = (Mach^-0.5)*Vvec(ii,1)*holder2(1) + ...
                   (Mach^-0.5)*Vvec(ii,2)*holder2(2) + ...
                   (Mach^-0.5)*Vvec(ii,3)*holder2(3);
%    Machuncer(ii) = sqrt(sum(holder2));

end

return