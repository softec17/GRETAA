function Machuncer = calculateMachUncerVinhEOM(Machlist,Pvec,statelist,gamma,soundspeedHist)
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
   Mach = Machlist(ii);
   state = statelist(ii,:);
   speedSound = soundspeedHist(ii);
%    speedSound = sqrt(gamma*state(11)/state(12));
   holder = makemat(Pvec(ii,:));
%    holder = holder.^0.5;
%    if isreal(holder)~=1
%        keyboard
%    end
   holder2 = holder(4,4).^2/speedSound;
%    + state(4)*sqrt(gamma/state(12))*0.5*holder(11,11).^2*state(11)^-0.5/speedSound^2 - ...
%        state(4)*sqrt(gamma*state(11))*(-0.5)*holder(12,12).^2*state(12)^-1.5/speedSound^2;
   realcheck = isreal(holder2);
   if realcheck ~= 1
       keyboard 
%        error(['Imaginary standard deviation found. ii = ',num2str(ii)])
   end
%    if isnan(holder2)==1
%        keyboard
%    end
%    holder2 = sqrt(real(holder2));     % Turn covariance into standard deviation
   Machuncer(ii) = holder2;
end

return