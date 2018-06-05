%% script to assess performance of group detection as tuning network from random to modular 
% Here in script 02: Communities detected by community detection algorithms (using Q >
%
%
% Change log:
% 5/6/2018: initial version
%
% Mark Humphries

clear all; close all;
addpath('../SyntheticModel/');
addpath('../Network_Spectra_Functions/');
addpath('../ZhangNewman2015/');

%% load networks from rejection script
fname = 'P_rejection_SyntheticEqual_NoNoise_20180605T100138';
load(['../Results/' fname])

%% clustering parameters
clusterpars.nreps = 20;
clusterpars.nLouvain = 1;

%% loop: make models and do rejection etc
n = sum(Model.N);
nBatch = size(Results.Time,2);

ResultsFields = {'normVIQmaxSpectra','normVIConsensusSpectra','normVILouvain','normVIMultiway','Time',...
                    'nGrpsLouvain','nGrpsMultiway','nGrpsQmaxSpectra','nGrpsConsensusSpectra'};
fieldsize = [nBatch,1]; % parfor cannot handle matrices of structs...
ClustResults = emptyStruct(ResultsFields,fieldsize);

blnP = autoParallel;  % set-up parallel processing with scaled number of cores.

parfor iB = 1:nBatch
% for iB = 1:nBatch
    
    for iP = 1:numel(Model.P_of_within)
        disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_within)) '; batch: ' num2str(iB)])
    
        tic 
        % Louvain on synthetic network
        [LouvCluster,LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Network(iP,iB).W,1,1);  % return 1st level of hierarchy only
        ClustResults(iB).normVILouvain(iP) = VIpartitions(LouvCluster{1}{1},group_membership) ./ log(numel(group_membership));
        ClustResults(iB).nGrpsLouvain(iP) = max(LouvCluster{1}{1});

        % multi-way spectra on synthetic network
        [bestPartition] = multiwaySpectCommDet(Network(iP,iB).W);
        ClustResults(iB).normVIMultiway(iP) = VIpartitions(bestPartition,group_membership) / log(numel(group_membership));
        ClustResults(iB).nGrpsMultiway(iP) = max(bestPartition);

         %% specified clustering on synthetic network
        if Results.SpectraWCMGroups(iP,iB) > 1
            [QmaxCluster,Qmax,ConsCluster,ConsQ,~] = ...
                                                        ConsensusCommunityDetect(Network(iP,iB).W,Network(iP,iB).ExpW,Results.SpectraWCMGroups(iP,iB),Results.SpectraWCMGroups(iP,iB),clusterpars.nreps);
            % quality of estimation of retained communities
            ClustResults(iB).normVIQmaxSpectra(iP)=VIpartitions(QmaxCluster,group_membership) / log(numel(group_membership));
            ClustResults(iB).normVIConsensusSpectra(iP)=VIpartitions(ConsCluster,group_membership) / log(numel(group_membership));
            ClustResults(iB).nGrpsQmaxSpectra(iP) = max(QmaxCluster);
            ClustResults(iB).nGrpsConsensusSpectra(iP) = max(ConsCluster);
        else
            QmaxCluster = []; Qmax = 0; ConsCluster = []; ConsQ = 0;
            ClustResults(iB).normVIQmaxSpectra(iP)=0;
            ClustResults(iB).normVIConsensusSpectra(iP)=0;
            ClustResults(iB).nGrpsQmaxSpectra(iP) = 0;
            ClustResults(iB).nGrpsConsensusSpectra(iP) = 0;
        end
        ClustResults(iB).Time(iP) = toc;
    end
end

%% save
save(['../Results/Clustering' fname],'ClustResults','clusterpars')


