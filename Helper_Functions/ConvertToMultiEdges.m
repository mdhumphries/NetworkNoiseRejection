function [A,conversion] = ConvertToMultiEdges(W,conversion)

% CONVERTTOMULTIEDGES convert weighted network to multi-edges
% [A,C_FINAL] = CONVERTTOMULTIEDGES(W,C) converts the weighted network in matrix W to
% a multi-edge network, given the conversion factor C. 
% The conversion is A_ij = round(W_ij*C), with C >= 1, turning any weights into integer
% counts of multi-edges.
%
% Set C = 'all' to automatically set smallest non-zero weight in the
% network to 1 (warning: this may create a large conversion factor, and so
% a network A with many multi-edges, slowing the null model generation).
%
% Will give warnings if the conversion factor will turn non-zero entries in
% W into zero entries in A (i.e. will delete links)
%
% Returns: 
% A, the matrix of multi-edges, same size as W.
% C_FINAL, the final value of C used for the network
%
% 7/9/2020: initial version, from existing code in POISSONSPARSEWCM
% Mark Humphries

minW = min(min(W(W>0)));  % minimum value of weight in the network
if strcmp(conversion,'all')
    % the scale so that minimum non-zero weight is 1
    conversion = 1./minW;
elseif ischar(conversion)       % conversion is a string, but not 'all'
    error('ConvertToMultiEdges:parameter','Unknown conversion option string specified')
elseif conversion < 0 || rem(conversion,1) ~= 0
    error('ConvertToMultiEdges:parameter','Conversion factor must be a positive integer')
end

% check if weights are already integers, then overwrite all options
if ~any(rem(W(:),1))  % then is integers for all weights
    conversion = 1;
end

% convert network into multi-edge version
A = W*conversion; 

% check how many entries will be set to 0 by conversion
idxs = find(A > 0 & A < 0.5);
nZeros = numel(idxs);
if nZeros > 0 warning('ConvertToMultiEdges:parameter',['Using a conversion factor of ' num2str(conversion) ' will set ' num2str(nZeros) ' existing links in your data network to zero']); end

% create integer edge counts
A = round(A);     



