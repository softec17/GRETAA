function b = specialInterp(t,a,time_a)
% Use linear interpolation algebraically to save comp. time
% a is a column matrix (i.e. number of rows greater than number of cols)
% time_a is a column vector with the same number of rows as a

length_a = max(size(a));

index = find(time_a<=t,1,'last');

if time_a(index) == t
    b = a(index,:);
elseif index == 0 || index == length_a
    b = interp1(time_a,a,t);
else
    slope = (a(index+1,:)-a(index,:))./(time_a(index+1)-time_a(index));
    b = a(index,:)+slope.*(t-time_a(index));
end

return