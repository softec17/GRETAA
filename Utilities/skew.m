function skewmat = skew(vec)
% Create a skew matrix used in cross products
%
% Input
% vec = vector of 3 x 1
%
% Output
% skewmat = skew matrix

skewmat = zeros(3);
skewmat(1,2) = -vec(3);
skewmat(1,3) = vec(2);
skewmat(2,3) = -vec(1);
skewmat(2,1) = vec(3);
skewmat(3,1) = -vec(2);
skewmat(3,2) = vec(1);

return