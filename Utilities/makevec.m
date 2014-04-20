% File: makevec
% By: Soumyo Dutta
% Modified: September 2009
% File takes a symmetric matrix and creates a vector for that
% Primarily used to make covariance matrices into vectors and vice versa

function vec = makevec(mat)
% Takes the upper triangle of a symmetric matrix and makes it into a vector
%
% Inputs
% mat = symmetric matrix
%
% Outputs
% vec = vector containing the upper triangle elements of a symmetric matrix

len_state = max(size(mat));

% Computes the index location in the covariance vector of the diagonal element 
% and end of row elements in the covariance matrix
diagonal = zeros(len_state,1);
rowends  = diagonal;
diagind = len_state:-1:1;
vec = zeros(sum(diagind),1);
diagonal(1) = 1;
rowends(1)  = len_state;
for ii=2:1:len_state
    diagonal(ii) = diagonal(ii-1)+diagind(ii-1);
    rowends(ii) = diagonal(ii)+(diagind(ii)-1);
end

% Convert the upper triangle of a covariance matrix to a vector
for ii = 1:1:len_state
    vec(diagonal(ii):rowends(ii)) = mat(ii,ii:end)';
end


return