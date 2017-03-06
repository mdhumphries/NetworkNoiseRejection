function Q = computeQ(C,B,m)

% COMPUTEQ modularity score of a multiple-partition 
% Q = COMPUTEQ(C,B,M) computes the modularity score Q, given:
%   C: a column vector of group membership IDs (integer values, one per
%   node)
%   B: the modularity matrix
%   M: the number of unique edges or sum of unique weights (i.e. that undirected networks have had each link counted once only)
%
%   Reference:  Newman, M. E. J. (2006) "Finding community structure in
%   networks using the eigenvectors of matrices". Phys Rev E, 74, 036104.
%
% Mark Humphries 2/3/2017

n = size(B,1);
grpIDs = unique(C);
ngrps = numel(grpIDs);

% construct S matrix of group membership: each column is a group
% See: Newman (2006) Eq 31            
S = zeros(n,ngrps);

for loop = 1:ngrps
    S(:,loop) = (C == grpIDs(loop));
end

% compute modularity
Q = trace(S' * B * S) / (2*m);  % note: assumes here that m is the sum of unique weights (i.e. that undirected networks have had each link counted once only)
