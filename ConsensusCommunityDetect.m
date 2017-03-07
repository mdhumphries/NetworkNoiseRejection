function [grps,Qmax,grpscon,Qcon,ctr,varargout] = ConsensusCommunityDetect(W,P,M,varargin)
 
% CONSENSUSCOMMUNITYDETECT partition signal network using eigenvectors of signal modularity matrix (with consensus)
%   [C,Qmax,Ccon,Qc,N,Q] = CONSENSUSCOMMUNITYDETECT(W,P,M) splits the
%   vertices of the nxn weighted signal network W into multiple groups, given
%   expected null model P, and maximum number of groups M.
%   
%   Returns: 
%       C: column vector indicating group membership for the parition with maximum
%           modularity
%       Qmax: the corresponding modularity score. 
%       Ccon: column vector indicating group membership for the consensus
%           clustering partition [see Notes]; 
%       Qcon: the corresponding is the modularity score for the consensus partition.
%       N:  the number of iterations until consensus was reached. 
%   
%   ...= CONSENSUSCOMMUNITYDETECT(...,N,DIMS) sets k-means to run  N times for each 
%       specified group size (default is 50); uses either 'all' (default)
%       or 'scaled' embedding dimensions for each tested K. See KMEANSSWEEP
%
%   [...,CLU] = CONSENSUSCOMMUNITYDETECT(...) where CLU is an optional output argument, 
%   returns every single clustering of the adjacency matrix A in the first
%   pass (i.e. before the consensus) - this is useful for further
%   post-processsing.
%
%   Notes: 
%   (0) Adjacency matrix: this can be weighted and directed. Directed networks are symmetrised as (W+W')/2 
%   When analysing time-series, this can be the similarity matrix. However, when
%   starting from a similarity matrix, ensure: 
%       (i) no self-loops - diagonal of A is all zeros; 
%       (ii) it's a similiarity matrix, not a correlation matrix: no
%       negative values
%   Warnings for both of these will be given
%
%   (1) This is a one-step multiple partition method, following up a
%   suggestion in Newman (2006) that all eigenvectors corresponding to positive
%   eigenvalues of the modularity matrix contain information about group
%   structure. The algorithm implemented takes the C such eigenvectors, and
%   uses k-means clustering on those eigenvectors to cluster the nodes into k = C+1 groups. 
%   A value for Q is computed for each k-means clustering (using the defined distance metrics).
%
%   Q = 1/2m * sum_ij(A_ij-P_ij) I(n_i,n_j) [where I(n_i,n_j) = 1 if nodes i,j
%   are in the same group, and 0 otherwise]
%   
%   (2) This is repeated for each C in 1:M, where M is number of positive
%   eigenvalues
%
%   (3) Consensus: this attempts to extract a stable set of groups that are robust to repeats
%   of the clustering process. All clusterings with Q>0 across all k-means variants and numbers of groups are
%   pooled. A consensus matrix is computed (Lancichinetti & Fortunato 2012): entry p_ij gives the
%   proportion of clusterings that placed nodes i and j in the same group.
%   The consensus matrix is then run through the one-step multiple parition
%   method (eigenvectors and k-means clustering). A new consensus matrix is
%   created. This is repeated until the distribution of p_ij has become
%   sufficiently bimodal, indicating that the groupings are stable. 

%   (5) For the community-detection algorithm, kmeans centres are initialised using the kmeans++ algorithm (Arthur & Vassilvitskii, 2007)
%
%   (6) Detection of bimodal distribution of consensus matrix now sets the
%   k-means initial centres dynamically to the 5th and 95th percentile of
%   the entries on the consensus matrix. 
%
%   References: 
%   (1) Newman, M. E. J. (2006) "Finding community structure in
%   networks using the eigenvectors of matrices". Phys Rev E, 74, 036104.
%
%   (2) Reichardt & Bornhaldt (2006) "Statistical mechanics of community detection".
%   Phys Rev E. 74, 016110
%   
%   (3) Lancichinetti, A. & Fortunato, S. (2012) Consensus clustering in complex networks.
%   Scientific Reports, 2, 336
%   
%   (4) Arthur, D. & Vassilvitskii, S. (2007) k-means++: the advantages of careful seeding. 
%   SODA '07: Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms, Society for Industrial and Applied Mathematics, 1027-1035
%
%   Mark Humphries 2/3/2017
nreps = 50;     % of each distance metric
dims = 'all';   % use all embedding dimensions for each k-means clustering

%% check if the passed matrix is a graph: catch common errors when passing a similarity matrix

% (1) no self-loops allowed
if ~all(diag(W)==0) 
    warning('Results likely unreliable: adjacency matrix has self-loops. Set diagonal to zero if no self-loops are needed.')
end

% (2) no negative links
x = sum(sum(W < 0));
if x > 0
    warning('Results likely unreliable: adjacency matrix has negative values')
end

%% set up options
if nargin >= 4
    if ~isempty(varargin{1}) nreps = varargin{1}; end
    if ~isempty(varargin{2}) 
        dims = varargin{2}; 
    end    
end

% % set up saving of each iteration
% blnSave = 0; % internal flag for setting saving of data
% if blnSave  
%     fname = ['Consensus_Iterations_' datestr(now,30)];
%     save(fname,'dims','nreps');  % save initial data to allow -append to work below
% end

%% internal parameters

nIDs = size(W,1);     % number of nodes of the weight matrix
m = sum(sum(W))/2;    % number of unique links (or total unique weights)

blnConverged = 0;       % stopping flag
ctr = 1;                % iterations of consensus

%% cluster signal network
B = W - P;          % initial modularity matrix, given data matrix W and specified null model P
[V,egs] = eig(B,'vector');
[~,ix] = sort(egs,'descend');    % sort into descending order
V = V(:,ix);                       % ditto the eigenvectors 

% C = kmeansSweep(V(:,1:M-1),2,M,nreps,dims);  % find groups in embedding dimensions: sweep from 2 to M
C = kmeansSweep(V(:,1:M),2,M,nreps,dims);  % find groups in embedding dimensions: sweep from 2 to M

for iQ = 1:size(C,2)
    Q(iQ) = computeQ(C(:,iQ),B,m); % compute modularity Q for each clustering
end

%% check for exit or store Qmax results
if isempty(C) || all(Q(:) <= 0)
    % then no groups detected; return empty
    grps = zeros(nIDs,1);
    Qmax = 0; 
    grpscon = zeros(nIDs,1);
    Qcon = 0;
    ctr = 0; 
    varargout{1} = C;
    blnConverged = 1;  % no answer
    return
else
    Qmax = max(Q);
    ix = find(Qmax == Q);
    grps = C(:,ix(1));  
    varargout{1} = C;
end
%% loop: consensus matrix and its clustering until converged

while ~blnConverged
    % make consensus
    Allowed = (Q > 0);       % only take consensus using Q-positive clusterings...
    CCons = makeConsensusMatrix(C(:,Allowed));
    
    %% check convergence
    [blnConverged,grpscon] = CheckConvergenceConsensus(CCons);
    
    %% sort out what to do next
    if blnConverged
        Qcon = computeQ(grpscon,B,m);  % compute Q, and exit
    else
        ctr = ctr+1; % increment consensus iteration counter
        if ctr > 50
            % do escape if not converging
            warning('Did not converge in 50 iterations - exiting without consensus answer')
            grpscon = [];
            Qcon = 0;
            blnConverged = 1;
        else
%             % find upper limit of groups - replicate this code when using
%             null model for consensus
%             [D,~,Mcons] = EmbedConsensusWishart(CCons);
%              % do k-means sweep using found M
%             C = kmeansSweep(D,2,Mcons,nreps,dims);  % find groups in embedding dimensions
           
            % do Laplacian on consensus matrix, using original M
            D = ProjectLaplacian(CCons,M);
            C = kmeansSweep(D,M,M,nreps,dims);  % find groups in embedding dimensions
             
            % compute Q
            Q = zeros(size(C,2),1);
            for iQ = 1:size(C,2)
                Q(iQ) = computeQ(C(:,iQ),B,m); % compute modularity Q for each clustering using original modularity matrix
            end
        end
    end
end



    
     




  







