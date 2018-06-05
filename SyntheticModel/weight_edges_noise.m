function [W,varargout] = weight_edges(A,N,S,alpha,P)

% WEIGHT_EDGES assign weights to edges given strength sequence
% W = WEIGHT_EDGES(A,N,S,ALPHA) creates a weight matrix for a network, given 
% a n*n adjacency matrix A (entries {0,1}), the set
% of sizes of each of the G modules N = [N1,N2,..,NG], and a n-length strength sequence
% S. Parameter ALPHA tunes the proportion of the node's strength that is
% assigned within and betweeen modules:
%       ALPHA = -1: all weights assigned between modules
%       ALPHA = 0: the proportion of weights follows the proportions of
%       links
%       ALPHA = 1: all weights assigned within module
%
% The function also needs:
%       P : the struct if probabilities used to make A (see WIRE_EDGES)
%
% Returns the n*n weight matrix W 
%
%
% Notes:
%   Computes the expected strength over all nodes as the mean of array S,
%   as this array can be generated from a variety of underlying
%   distributions. 
%
% Change log:
% 07/08/2017: original version - uses ALPHA to set S(in) = ALPHA  * S;
% 22/05/2018: ALPHA tunes the strength distribution according to the link
% distribution
%
% Mark Humphries

% check alpha
if alpha < -1 || alpha > 1
    error('weight_edges:parameter','ALPHA is out of range')
end

% check whether there are strengths to place...
ixpairs = find((triu(A,0)>0));     % linear indices of linked pairs for upper triangular matrix
nLinks = round(sum(S/2)) - numel(ixpairs);  % number of links left to place: total links - [already placed]

if nLinks < 0
    warning('Sampled strengths do not exceed existing degree distribution')
    W = A;
    return
end

%% calculate baseline proportion of links within and between groups
n = sum(N);                         % number of nodes in network
E_in_group = (N-1)/2 .* P.in;       % expected number of unique links in group per node
E_out_group = (n-N)/2 .* P.between; % expected number of unique links outside of group

Sratio = E_in_group ./ (E_in_group + E_out_group);  % proportion of links within groups

%% define tuning of strengths
max_alpha = 1./Sratio - 1;   % to give S(in) == S [all strengths within modules]

if alpha > 0 
    alpha_star = alpha .* max_alpha;
else
    alpha_star = zeros(1,numel(N)) + alpha;
end


%% group membership indices, and alpha per node...
G = zeros(n,1);

% group membership
Nsum = [0 cumsum(N)];
G = zeros(n,1); alpha_per_node = zeros(n,1); Sratio_per_node = zeros(n,1);
for iG = 1:numel(Nsum)-1
    G(Nsum(iG)+1:Nsum(iG+1)) = iG;
    alpha_per_node(Nsum(iG)+1:Nsum(iG+1)) = alpha_star(iG);
    Sratio_per_node(Nsum(iG)+1:Nsum(iG+1)) = Sratio(iG);
end

Grow = repmat(G',n,1);  % assign group IDs per row    
Gcol = repmat(G,1,n);   % assign group IDs per column

        
%% do strength distribution in two passes
% note that this truncates extremely high strengths, as we are giving
% Poisson distribution in two blocks
blnWithin = Grow == Gcol;       % boolean masks for within modules
blnBetween = Grow~=Gcol;        % boolean mask for between modules

% pass 1: within modules
S_within = (1+alpha_per_node) .* Sratio_per_node .* S; 
tempW = MakeWeights(A,S_within,blnWithin);

% pass 2: between modules
S_between = S - S_within;
tempWb = MakeWeights(A,S_between,blnBetween);

W = tempW + tempWb;  % combine the two weight matrices

varargout{1} = blnWithin;
varargout{2} = blnBetween;

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

