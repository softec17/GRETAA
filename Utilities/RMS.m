function RMSVal = RMS(x)

n_x = length(x);
RMSVal = sqrt(sum(x.^2)./n_x);

return