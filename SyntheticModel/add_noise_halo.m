function W = add_noise_halo(Wcore,Pnoise,Snoise)

% ADD_NOISE_HALO add halo of noise nodes to core network
% W = ADD_NOISE_HALO(C,P,S) adds a halo of noise nodes to a weighted block 
% model network, given:
%   C: the nxn weight matrix of the block model (the "core")
%   P: the probability of links from noise nodes
%   S: h-length sampled strength sequence for the noise nodes
%
% Returns the (n+h) x (n+h) weight matrix W, with the last h rows and 
% columns defining the noise halo
%
%
% Change log:
% 05/06/2018: original version
%
% Mark Humphries

n = size(Wcore,1);   % nodes in core
h = length(Snoise);  % nodes in halo

Grow = repmat(1:n+h,n+h,1);   % assign node IDs per row    
Gcol = repmat(1:n+h,1,n+h);   % assign node IDs per column
blnCoreNoise = Grow > n && Gcol <= n;       % boolean mask for links between Core and Halo
blnNoiseCore = Grow <= n && Gcol > n;       % boolean mask for links between Halo and Core
blnNoiseNoise = Grow > n && Gcol > n;       % boolean mask for links within Halo

blnMask = blnCoreNoise + blnNoiseCore + blnNoiseNoise;

%% add edges
Plink = blnMask.*Pnoise;  % matrix of probabilities per link

A = real(rand(n,n) < triu(Plink,1));    % randomly chose edges in upper triangle
A = A + A'; % undirected

%% weight edges
tempW = zeros(numel(Snoise));
ixpairs = find((triu(A,0)>0) .* blnMask);     % linear indices of linked pairs for upper triangular matrix
tempW(ixpairs) = 1;  % entries from A
nLinks = round(sum(Snoise/2)) - numel(ixpairs);  % number of links left to place: total links - [already placed]

if nLinks > 0 % then there are edges left to place     
    % get as (i,j)
    [irow,jcol] = ind2sub([n,n],ixpairs);
    % get P(link): 
    Plink = Snoise(irow) .* Snoise(jcol);
    Plink = Plink ./ sum(Plink); % P(link is placed between each pair)

    lambda = nLinks .* Plink; % expected number of links

    Nlink = poissrnd(lambda);  % Poisson random number of links made

    tempW(ixpairs) = tempW(ixpairs) + Nlink; % add to existing links
end

tempW = tempW + tempW';  % make symmetric


%% add together
W = Wcore + Whalo;

