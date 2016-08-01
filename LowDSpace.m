function  [Dspace,ixpos,Dn,mci] = LowDSpace(B,Emodel,I)

% LOWDSPACE find and return low-dimensional axes for network
% [D,X,N,M] = LOWDSPACE(B,E,I) finds the P low-dimensional axes for the network, given: 
%       B: the (nxn) modularity matrix of the data network, defined using a null model (e.g Weighted Configuration Model)
%       E: the null-model eigenvalue distribution (n x #repeats of model) (from e.g. WeightedConfigModel) 
%       I: specified confidence interval on the maximum eigenvalue
%       (scalar)[enter 0 to just use the mean estimate]
%
%  Returns:
%       D: the n x P matrix of retained eigenvectors
%       X: the corresponding set of indices into the full eigenvector matrix
%       N: the number of eigenvectors retained
%       M: [m CI] mean and confidence interval on the maxium eigenvalue
%
%   Notes: a helper function for NodeRejection, but also useful for calling
%   in its own right to just obtain the low-D space
%
%  ChangeLog:   
%   28/7/2016: initial version
%    1/8/2016: input check; implemented maximum eigenvalue rejection
%
% Mark Humphries

n = size(B,1);
if size(Emodel,1) ~= n
    error('Eigenvalue matrix should be n (nodes) x N (repeats)')
end

[V,egs] = eig(B,'vector');  % eigenspectra of data modularity matrix
[egs,ix] = sort(egs,'descend'); % sort eigenvalues into descending order 
V = V(:,ix);  % sort eigenvectors accordingly

% option 1: compute mean and CI over largest eigenvalue
mx = max(Emodel); 
M = mean(mx);
CIs = CIfromSEM(std(mx),size(Emodel,2),I);
bnds = M+CIs; % confidence interval on maximum eigenvalue for null model

% option 2: just compute pooled distribution, and return bounds
% prctI = [I/2*100 100-I/2*100]; % rejection interval as symmetric percentile bounds
% bnds = max(prctile(Emodel(:),prctI)); % confidence interval on eigenvalue distribution for null model

% return dimensions
ixpos = find(egs >= bnds); % eigenvalues exceeding these bounds

% return answers
Dn = numel(ixpos);   % number of retained dimensions          
Dspace = V(:,ixpos);  % axes of retained dimensions
mci = [M CIs];

