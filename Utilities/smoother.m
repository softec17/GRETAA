% Smoother.m
% Smoothing algorithm
% Takes the output from the forward pass and the backward pass of the
% Extended Kalman Filter and combines the result for a best estimate

% By: Soumyo Dutta
% Modified: September 2009

disp('Running smoothing procedure')
fwd_x = fwdstate.state;
fwd_Pvec = fwdstate.Pvec;
bck_x = bckstate.state;
bck_Pvec = bckstate.Pvec;
% Time vector: should be same for both bck and fwd 
% except bck will have the time vector in descending order
smthstate.time = fwdstate.time;
time_smooth = smthstate.time;
len_time = length(time_smooth);
len_state = min(size(fwd_x));
len_Pvec = len_state*(len_state+1)/2;

% Preallocate vectors for speed
smthstate.state = zeros(len_time,len_state);
smthstate.Pvec = zeros(len_time,len_Pvec);

%% Calculate smoothed results
j_s = [len_time:-1:1];
for i = 1:1:len_time
   Pf = makemat(fwd_Pvec(i,:));             %diag(fwd_std(:,i).^2);
%    Pb = makemat(bck_Pvec(len_time-i+1,:));  %diag(bck_std(:,len_time-i+1).^2); 
   Pb = makemat(bck_Pvec(j_s(i),:));  %diag(bck_std(:,len_time-i+1).^2); 
   add_P = inv(Pf)+inv(Pb);
   P_smooth = inv(add_P);
%    keyboard
%    x_smooth = add_P\(Pf\fwd_x(i,:)'+Pb\bck_x(len_time-i+1,:)');
%    x_smooth = add_P\(Pf\fwd_x(i,:)'+Pb\bck_x(j_s(i),:)');
   x_smooth = P_smooth\(Pf\fwd_x(i,:)' + Pb\bck_x(j_s(i),:)');

   smthstate.state(i,:) = x_smooth';
   smthstate.Pvec(i,:) = makevec(P_smooth)';
end
