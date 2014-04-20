function normvec = norm_vec(vec)
% Function takes the norm of the components of a vector at every timeline
%
% Inputs
% vec = vector with components
%
% Outputs
% normvec = norm of the components

[r,c] = size(vec);
if r<c
    normvec = zeros(c,1);
    for ii = 1:1:c
        normvec(ii) = norm(vec(:,ii));
    end
else
    normvec = zeros(r,1);
    for ii = 1:1:r
        normvec(ii) = norm(vec(ii,:));
    end
end

return