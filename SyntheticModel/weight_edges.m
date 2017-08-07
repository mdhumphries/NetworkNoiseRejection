function W = weight_edges(A,N,S,alpha)

% WEIGHT_EDGES assign weights to edges given strength sequence
% W = WEIGHT_EDGES(A,N,S,ALPHA) creates a weight matrix for a network, given 
% a n*n adjacency matrix A (entries {0,1}), the set
% of sizes of each of the G modules N = [N1,N2,..,NG], and a n-length strength sequence
% S. Parameter ALPHA controls the proportion of the node's strength that is
% assigned within and betweeen modules.
%
% Returns the n*n weight matrix W 
%
% Mark Humphries 7/8/2017

n = sum(N);                      % number of nodes

%% group membership indices
G = zeros(n,1);

% group membership
Nsum = [0 cumsum(N)];
G = zeros(n,1);
for iG = 1:numel(Nsum)-1
    G(Nsum(iG)+1:Nsum(iG+1)) = iG;
end

Grow = repmat(G',n,1);  % assign group IDs per row    
Gcol = repmat(G,1,n);   % assign group IDs per column
blnWithin = Grow == Gcol;
blnBetween = Grow~=Gcol;

%% do WCM model in two passes

% pass 1: within modules
tempS = S .* alpha;
tempW = MakeWeights(A,tempS,blnWithin);

% pass 2: between modules
tempS = S .* (1-alpha);
tempWb = MakeWeights(A,tempS,blnBetween);

W = tempW + tempWb;  % combine the two weight matrices

%% Poisson function
function tempW = MakeWeights(A,S,blnMask)
    % blnMask: n*n boolean matrix of included nodes
    n = numel(S);
    tempW = zeros(numel(S));
    ixpairs = find((triu(A,0)>0) .* blnMask);     % linear indices of linked pairs for upper triangular matrix
    tempW(ixpairs) = 1;  % entries from A
    nLinks = round(sum(S/2)) - numel(ixpairs);  % number of links left to place: total links - [already placed]

    if nLinks > 0 % then there are edges within-modules left to place     
        % get as (i,j)
        [irow,jcol] = ind2sub([n,n],ixpairs);
        % get P(link): 
        Plink = S(irow) .* S(jcol);
        Plink = Plink ./ sum(Plink); % P(link is placed between each pair)

        lambda = nLinks .* Plink; % expected number of links

        Nlink = poissrnd(lambda);  % Poisson random number of links made
        
        tempW(ixpairs) = tempW(ixpairs) + Nlink; % add to existing links
    else
        warning('Sampled strengths do not exceed existing degree distribution')
    end

    tempW = tempW + tempW';  % make symmetric

