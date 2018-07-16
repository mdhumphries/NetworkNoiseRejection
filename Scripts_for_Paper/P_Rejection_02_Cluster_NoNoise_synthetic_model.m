%% script to assess performance of group detection as tuning network from random to modular 
% Here in script 02: Communities detected by community detection algorithms (using Q >
%
%
% Change log:
% 5/6/2018: initial version
%
% Mark Humphries

clearvars Network Results ClustResults

addpath('../Network_Analysis_Functions/');
addpath('../ZhangNewman2015/');
addpath('../Helper_Functions/')

%% load networks from rejection script
% comment out filename to run batch script
% fname = 'P_rejection_SyntheticEqual_NoNoise_20180702T124726';  % full set of 100 networks per P(within) level
fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';

load([fpath fname])

%% clustering parameters
clusterpars.nreps = 50;
clusterpars.nLouvain = 1;

%% loop: make models and do rejection etc
n = sum(Model.N);
nBatch = size(Results.Time,2);

ResultsFields = {'normVIQmaxFull','normVIConsensusFull','normVIQmaxSparse','normVIConsensusSparse','normVILouvain','normVIMultiway','normVIMultiwayQ','Time',...
                    'nGrpsLouvain','nGrpsMultiway','nGrpsMultiwayQ','nGrpsQmaxFull','nGrpsConsensusFull','nGrpsQmaxSparse','nGrpsConsensusSparse'};
fieldsize = [nBatch,1]; % parfor cannot handle matrices of structs...
ClustResults = emptyStruct(ResultsFields,fieldsize);
LoopResults = emptyStruct(ResultsFields,fieldsize);

blnP = autoParallel;  % set-up parallel processing with scaled number of cores.

parfor iB = 1:nBatch
% for iB = 1:nBatch
    
    for iP = 1:numel(Model.P_of_within)
        disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_within)) '; batch: ' num2str(iB)])
    
        tic 
        % Louvain on synthetic network
        [LouvCluster,LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Network(iP,iB).W,1,1);  % return 1st level of hierarchy only
        [~,ClustResults(iB).normVILouvain(iP)] = VIpartitions(LouvCluster{1}{1},group_membership);
        ClustResults(iB).nGrpsLouvain(iP) = max(LouvCluster{1}{1});

        % multi-way spectra on synthetic network
        [bestPartition,maxQPartition] = multiwaySpectCommDet(Network(iP,iB).W,Results.SpectraSparseWCM.Groups(iP,iF,iB)*2);
        [~,ClustResults(iB).normVIMultiway(iP)] = VIpartitions(bestPartition,group_membership);
        ClustResults(iB).nGrpsMultiway(iP) = max(bestPartition);
        [~,ClustResults(iB).normVIMultiwayQ(iP)] = VIpartitions(maxQPartition,group_membership);
        ClustResults(iB).nGrpsMultiwayQ(iP) = max(maxQPartition);
           
        %% specified clustering on synthetic network
        % (1) sparse WCM
        LoopResults = clusterLowDNetwork(Network(iP,iB).W,Network(iP,iB).ExpW,Results.SpectraSparseWCM.Groups(iP,iB),Results.SpectraSparseWCM.Groups(iP,iB),clusterpars.nreps,group_membership);
        [Clust(iB).normVIQmaxFull(iP),Clust(iB).normVIConsensusFull(iP),Clust(iB).nGrpsQmaxFull(iP),Clust(iB).nGrpsConsensusFull(iP)] ...
            = deal(LoopResults.normVIQmaxSpectra,LoopResults.normVIConsensusSpectra,LoopResults.nGrpsQmaxSpectra,LoopResults.nGrpsConsensusSpectra);

        % (2) full WCM
        P = expectedA(Network(iP,iB).W);  % construct P on the fly, as is quick...
        LoopResults = clusterLowDNetwork(Network(iP,iB).W,P,Results.SpectraFullWCM.Groups(iP,iB),Results.SpectraFullWCM.Groups(iP,iB),clusterpars.nreps,group_membership);
        [Clust(iB).normVIQmaxFull(iP),Clust(iB).normVIConsensusFull(iP),Clust(iB).nGrpsQmaxFull(iP),Clust(iB).nGrpsConsensusFull(iP)] ...
            = deal(LoopResults.normVIQmaxSpectra,LoopResults.normVIConsensusSpectra,LoopResults.nGrpsQmaxSpectra,LoopResults.nGrpsConsensusSpectra);

            
        ClustResults(iB).Time(iP) = toc;
    end
end

%% save
save([fpath 'Clustering' fname],'ClustResults','clusterpars')


