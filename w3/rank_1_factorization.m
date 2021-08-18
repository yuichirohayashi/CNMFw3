function [spmat, tempmat] = rank_1_factorization(Y, iter)

% Rank-1 factorization
% approximate Y as spmat*tempmat

if nargin < 2
    iter=1;
end

spmat = ones(size(Y,1),1)/sqrt(size(Y,1));
for iter=1:iter
    tempmat = spmat'*Y;
    tempmat = tempmat/norm(tempmat(:));
    spmat = Y*tempmat';
    if iter<iter
        spmat = spmat/norm(spmat(:));
    end
end
