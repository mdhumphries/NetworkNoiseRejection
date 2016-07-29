function  [Dspace,ixpos,Dn] = LowDSpace(B,Emodel,I)

% LOWDSPACE find and return low-dimensional axes for network
% [D,X,N] = LOWDSPACE(B,E,I) finds the P low-dimensional axes for the network, given: 
%       B: the (nxn) modularity matrix of the data network, defined using a null model (e.g Weighted Configuration Model)
%       E: the null-model eigenvalue distribution (n x #repeats of model) (from e.g. WeightedConfigModel) 
%       I: specified rejection interval (scalar)
%
%  Returns:
%       D: the n x P matrix of retained eigenvectors
%       X: the corresponding set of indices into the full eigenvector matrix
%       N: the number of eigenvectors retained
%   
%   Notes: a helper function for NodeRejection, but also useful for calling
%   in its own right to just obtain the low-D space
%
%  ChangeLog:   
%   28/7/2016: initial version
%
% Mark Humphries

[V,egs] = eig(B,'vector');  % eigenspectra of data modularity matrix
[egs,ix] = sort(egs,'descend'); % sort eigenvalues into descending order 
V = V(:,ix);  % sort eigenvectors accordingly

% rejection interval on distribution
prctI = [I/2*100 100-I/2*100]; % rejection interval as symmetric percentile bounds
bnds = prctile(Emodel(:),prctI); % confidence interval on eigenvalue distribution for null model
ixpos = find(egs >= bnds(2)); % eigenvalues exceeding these bounds

% return answers
Dn = numel(ixpos);   % number of retained dimensions          
Dspace = V(:,ixpos);  % axes of retained dimensions

