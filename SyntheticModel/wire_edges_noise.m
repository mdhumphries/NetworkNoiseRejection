function A = wire_edges_noise(N,F,P)

% WIRE_EDGES_NOISE create planted partition block model for edges (undirected)
% A = WIRE_EDGES_NOISE(N,F,P)creates a planted partition block model in a 
% halo of noise nodes, given the set of sizes of each of the G modules 
% N = [N1,N2,..,NG], the fraction of additional nodes that make up the noise halo,
% (F=1 is half all nodes are noise), and the probabilities of edges between them in struct P:
%       P.in : probability of within-module edges
%       P.between: probability of between-module edges
%       P.noise : probability of edges from/to/between noise nodes
%
% Returns the n*n adjacency matrix A (where n = sum(N))
%
% Change log:
% 06/06/2018: initial version  
%
% Mark Humphries

n = sum(N); 
n_noise = n * F;
T_nodes = n + n_noise;
G = zeros(T_nodes,1);  % halo coded as 0

% group membership - of core nodes
Nsum = [0 cumsum(N)];
for iG = 1:numel(Nsum)-1
    G(Nsum(iG)+1:Nsum(iG+1)) = iG;
end

% N_in = sum(N.^2) - n;  % number of within-module connections no self-connections
% N_out = n^2 - N_in;    % numner of between-module connections

Grow = repmat(G',T_nodes,1);  % assign group IDs per row    
Gcol = repmat(G,1,T_nodes);   % assign group IDs per column
Plink = (Grow == Gcol & Gcol > 0 & Grow > 0).*P.in + (Grow~=Gcol & Gcol > 0 & Grow > 0).*P.between + (Grow==0 | Gcol==0) .* P.noise;  % matrix of probabilities per link

A = real(rand(T_nodes,T_nodes) < triu(Plink,1));    % randomly chose edges in upper triangle
A = A + A'; % undirected