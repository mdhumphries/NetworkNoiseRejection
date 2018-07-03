% script to test VI on clusters + noise

addpath ../Helper_Functions/
addpath ../Network_Analysis_Functions/

% ground-truth parameters
N = [100 100 100];
fnoise = 0.5;

G = numel(N);
Tcluster = sum(N);
Tnoise = Tcluster * fnoise;
T = Tcluster + Tnoise;

% label ground-truth clusters
cumul_N = [0 cumsum(N)];
for iG = 1:G
    clusterIDs(1+cumul_N(iG):cumul_N(iG+1)) = iG;
end

%% 1. ground-truth is all noise as separate 1 node groups

noiseIDs = G+1:Tnoise+G;

partition = [clusterIDs'; noiseIDs(randperm(Tnoise))']; % perfect clusters, plus permuted one node groups

VI = VIpartitions([clusterIDs'; noiseIDs'],partition);  % should be 0!

%% 2. partition is random permutation of separate 1 node groups with clusters (of exactly correct size)

partition2 = partition(randperm(T));
VI2 = VIpartitions([clusterIDs'; noiseIDs'],partition2) / log2(T);

%% 3. ground-truth is all noise as own group
noiseIDs2 = ones(Tnoise,1) + G;

VIrandomnoise = VIpartitions([clusterIDs'; noiseIDs2],partition) / log2(T);  % noise as 1 group vs perfect-clusters  & one-per-group-noise

VIallrandom = VIpartitions([clusterIDs'; noiseIDs2],partition2) / log2(T); % noise as 1 group vs randomised-clusters+noise together

VInogroups = VIpartitions([clusterIDs'; noiseIDs2],randperm(T)) / log2(T);  % cluster+noise versus random 1-per-group
VInogroups = VIpartitions([clusterIDs'],randperm(Tcluster)) / log2(T);  % cluster-only versus random 1-per-group


%% 3. partition is missing nodes from ground truth....
% remove rejected nodes from ground truth too

