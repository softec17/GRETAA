% File: makemat
% By: Soumyo Dutta
% Modified: September 2009
% File takes a vector and makes a symmetric matrix for that vector
% Primarily used to make covariance matrices into vectors and vice versa

function mat = makemat(vec)
% Takes the vector representation of the upper triangle of a symmetric
% matrix and converts it into matrix form
% Inputs
% vec = vector containing the upper triangle elements of a symmetric matrix
%
% Outputs
% mat = symmetric matrix


% Note: N*(N+1)/2 = sum of consecutive numbers from N to 1
% Formula used to find length of vec; If vec_length known, solve for N
% N = number of states
len_vec = length(vec);
len_state = -0.5 + 0.5*sqrt(1+4*2*len_vec); 

% Computes the index location in the covariance vector of the diagonal element 
% and end of row elements in the covariance matrix
diagonal = zeros(len_state,1);
rowends  = diagonal;
diagind = len_state:-1:1;
diagonal(1) = 1;
rowends(1)  = len_state;
for ii=2:1:len_state
    diagonal(ii) = diagonal(ii-1)+diagind(ii-1);
    rowends(ii) = diagonal(ii)+(diagind(ii)-1);
end

mat = zeros(len_state);        % Preallocate space for speed
% Creates the upper-triangle of the eventually symmetric, covariance matrix
for ii = 1:1:len_state                  
    row = vec(diagonal(ii):rowends(ii))';
    mat(ii,ii:end) = row;
end
mat = mat + mat';    % Creates symmetric covariance based on 
                                    % the upper triangle
for ii = 1:1:len_state
    mat(ii,ii) = mat(ii,ii)/2;    
    % Process to make symmetric matrix double the diagonal elements. This
    % rectifies that
end

return