function [W] = weight_edges_noise(A,S)

% WEIGHT_EDGES_NOISE assign weights to edges given strength sequence
% W = WEIGHT_EDGES_NOISE(A,S) creates a weight matrix for a network, given 
% a n*n adjacency matrix A (entries {0,1}), and a n-length strength sequence
% S. 
%
% Returns the n*n weight matrix W 
%
%
% Notes:
%   1. The generative model distinguishes between in-module and out-module
%   connections; the out-module includes both between-module and noise halo
%   connections
%
% Change log:
% 06/06/2018: cloned original function to add handling of noise halo; 
%             no option for ALPHA here, as it is unclear how to define this 
%             parameter for the noise halo  
%
% Mark Humphries

% check whether there are strengths to place...
ixpairs = find((triu(A,0)>0));     % linear indices of linked pairs for upper triangular matrix
nLinks = round(sum(S/2)) - numel(ixpairs);  % number of links left to place: total links - [already placed]

if nLinks < 0
    warning('Sampled strengths do not exceed existing degree distribution')
    W = A;
    return
end
        
%% do strength distribution 
blnMask = triu(A,0)>0;       % boolean masks for within modules
W = MakeWeights(A,S,blnMask);


%% Poisson function
function tempW = MakeWeights(A,S,blnMask)
    % blnMask: n*n boolean matrix of included nodes
    n = numel(S);
    tempW = zeros(numel(S));
    ixpairs = find((triu(A,0)>0) .* blnMask);     % linear indices of linked pairs for upper triangular matrix
    tempW(ixpairs) = 1;  % entries from A
    nLinks = round(sum(S/2)) - numel(ixpairs);  % number of links left to place: total links - [already placed]
    
    
    if nLinks > 0 % then there are edges left to place     
        % get as (i,j)
        [irow,jcol] = ind2sub([n,n],ixpairs);
        % get P(link): 
        Plink = S(irow) .* S(jcol);
        Plink = Plink ./ sum(Plink); % P(link is placed between each pair)

        lambda = nLinks .* Plink; % expected number of links

        Nlink = poissrnd(lambda);  % Poisson random number of links made
        
        tempW(ixpairs) = tempW(ixpairs) + Nlink; % add to existing links
    end

    tempW = tempW + tempW';  % make symmetric

