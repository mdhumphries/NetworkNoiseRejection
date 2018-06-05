function A = wire_edges(N,P)

% WIRE_EDGES create planted partition block model for edges (undirected)
% A = WIRE_EDGES(N,P)creates a planted partition block model, given the set
% of sizes of each of the G modules N = [N1,N2,..,NG], and the
% probabilities of edges between them in struct P:
%       P.in : probability of within-module edges
%       P.between: probability of between-module edges
%
% Returns the n*n adjacency matrix A (where n = sum(N))
%
% Mark Humphries 3/8/2017

n = sum(N); 
G = zeros(n,1);

% group membership
Nsum = [0 cumsum(N)];
G = zeros(n,1);
for iG = 1:numel(Nsum)-1
    G(Nsum(iG)+1:Nsum(iG+1)) = iG;
end

% N_in = sum(N.^2) - n;  % number of within-module connections no self-connections
% N_out = n^2 - N_in;    % numner of between-module connections

Grow = repmat(G',n,1);  % assign group IDs per row    
Gcol = repmat(G,1,n);   % assign group IDs per column
Plink = (Grow == Gcol).*P.in + (Grow~=Gcol).*P.between;  % matrix of probabilities per link

A = real(rand(n,n) < triu(Plink,1));    % randomly chose edges in upper triangle
A = A + A'; % undirected