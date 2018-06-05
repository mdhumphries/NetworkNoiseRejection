%% script to assess performance on equal-sized groups 
% 
% Mark Humphries 7/8/2017

clear all; close all;
addpath('../SyntheticModel/');
addpath('../Network_Spectra_Functions/');
addpath('../ZhangNewman2015/');
addpath('../Helper_Functions/');

%% load synthetic model data
fname = 'SyntheticEqual_NoNoise_20180604T150903';
load(['../Results/' fname)

%% clustering parameters

clusterpars.nreps = 100;
clusterpars.nLouvain = 5;


%% loop: make models and do rejection etc
n = sum(Model.N);
nBatch = size(Results.Time,3);

FullFields = {'QmaxCluster','Qmax','ConsCluster','ConsQ','LouvCluster','LouvQ','VI_Louvain'};
ResultsFields = {'nVIFull_QmaxSpectra','nVIFull_ConsensusSpectra','nVIFull_LouvainMin','nVIFull_LouvainMax','nVIFull_Multiway','Time'};
fieldsize = [numel(Model.P_of_within),numel(Model.alpha_range),nBatch];
ClustResults = emptyStruct(ResultsFields,fieldsize);

% we can parallelise this
for iP = 1:numel(Model.P_of_within)
    for iA = 1:numel(Model.alpha_range)
        disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_within)) '; alpha ' num2str(iA) '/' num2str(numel(Model.alpha_range)) ])

        parfor iB = 1:nBatch
            tic   
            Full = emptyStruct(FullFields,[]);

            %% do all 3 detection methods on Signal?

            %% and all again, on full matrix, assigning noise to its own community
            if Data.Dn > 0
                [Full.QmaxCluster,Full.Qmax,Full.ConsCluster,Full.ConsQ,~] = ...
                                                            ConsensusCommunityDetect(Network(iP,iA,iB).A,Network(iP,iA,iB).ExpA,Results.SpectraWCMGroups(iP,iA,iB),Results.SpectraWCMGroups(iP,iA,iB));
            else
                Full.QmaxCluster = []; Full.Qmax = 0; Full.ConsCluster = []; Full.ConsQ = 0;
            end
            % quality of estimation of retained communities
            ClustResults(iP,iA,iB).nVIFull_QmaxSpectra=VIpartitions(Full.QmaxCluster,group_membership) / log(numel(group_membership));
            ClustResults(iP,iA,iB).nVIFull_ConsensusSpectra=VIpartitions(Full.ConsCluster,group_membership) / log(numel(group_membership));

            [Full.LouvCluster,Full.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Network(iP,iA,iB).A,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
            for j = 1:clusterpars.nLouvain
                CLou = Full.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
                Full.VI_Louvain(j) = VIpartitions(CLou,group_membership) ./ log(numel(group_membership));
            end
            ClustResults(iP,iA,iB).nVIFull_LouvainMin = min(Full.VI_Louvain);
            ClustResults(iP,iA,iB).nVIFull_LouvainMax = max(Full.VI_Louvain);

            % multi-way spectra
            [bestPartition] = multiwaySpectCommDet(Network(iP,iA,iB).A);
            ClustResults(iP,iA,iB).nVIFull_Multiway = VIpartitions(bestPartition,group_membership) / log(numel(group_membership));

            ClustResults(iP,iA,iB).Time = toc
        end
    end
end

%% save
save(['../Results/Clustering' fname],'ClustResults','clusterpars')


