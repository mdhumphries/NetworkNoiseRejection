% script to test Otsu vs k-means on pathological consensus matrices
% matrices from example Ca2+ imaging correlation matrix, from Peron data
% problem appears to be: matrix itself is converging, but the convergence
% detection (k-means) is not working.
%
% So: test Otsu

clear all; close all;

addpath ../Helper_Functions/
addpath ../Network_Analysis_Functions/

load CCons_10.mat

for i = 1:numel(CCons_10)
    
    [blnTrans(i),grpscon{i},theta(i)] = CheckConvergenceConsensus(CCons_10{i});
end

ixs = find(blnTrans);
n = size(CCons_10{1},1); 

for i = 1:numel(ixs)
    for j = i+1:numel(ixs)
        V(i,j) = VIpartitions(grpscon{ixs(i)},grpscon{ixs(j)}) ./ log(n);
    end
end

% test on actual network...

load ../Results//Rejected_Lesmis.mat
nreps = 200;

% construct new null model
P = Data.ExpA(Data.ixSignal_Final,Data.ixSignal_Final); % extract relevant part of null model
[Otsu.QmaxCluster,Otsu.Qmax,Otsu.ConsCluster,Otsu.ConsQ,ctr] = ...
                      ConsensusCommunityDetect(Data.Asignal_final,P,1+Data.Dn,1+Data.Dn,nreps);

load ../Results/Clustered_Lesmis.mat
ntest = size(P,1);
Vtest = VIpartitions(Otsu.ConsCluster,Connected.ConsCluster) ./ log(ntest)

% test again
load ../Results//Rejected_polblogs.mat
nreps = 200;

P = Data.ExpA(Data.ixSignal_Final,Data.ixSignal_Final); % extract relevant part of null model
[Otsu.QmaxCluster,Otsu.Qmax,Otsu.ConsCluster,Otsu.ConsQ,ctr] = ...
                      ConsensusCommunityDetect(Data.Asignal_final,P,1+Data.Dn,1+Data.Dn,nreps);

load ../Results/Clustered_polblogs.mat
ntest = size(P,1);
Vtest = VIpartitions(Otsu.ConsCluster,Connected.ConsCluster) ./ log(ntest)

