%% script to assess performance of group detection as tuning network from random to modular 
% Here in script 02: Communities detected by community detection algorithms
% on noise halo
%
%
% Change log:
% 08/06/2018: initial version
% 04/07/2018 : added comparison to ground-truths
% 16/07/2018 : updated to proper Full WCM
%
% Mark Humphries

clearvars Network Results ClustResults
addpath('../Network_Analysis_Functions/');
addpath('../ZhangNewman2015/');
addpath('../Helper_Functions/')

%% load networks from rejection script
% comment out filename to run batch script
% fname = 'P_rejection_SyntheticEqual_Noise_20180611T132723';  % full set of 100 networks per P(within) level
fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';

load([fpath fname])

%% clustering parameters
clusterpars.nreps = 50;
clusterpars.nLouvain = 1;
clusterpars.maxMultiway = 20;

%% define ground-truths
G = numel(Model.N);
Tcluster = sum(Model.N);

% ground truths
Nsum = [0 cumsum(Model.N)];
group_membership = zeros(Tcluster,1);
for iG = 1:numel(Nsum)-1
    group_membership(Nsum(iG)+1:Nsum(iG+1)) = iG;
end 

for iF = 1:numel(Model.F_noise)
    Tnoise = Tcluster * Model.F_noise(iF);

    % (1) noise nodes in each group by themselves
    noiseIDs = G+1:Tnoise+G;
    Partition(iF).owngroups = [group_membership; noiseIDs'];

    % (2) noise nodes in their own groups
    noiseIDs = ones(Tnoise,1) + G;
    Partition(iF).onegroup = [group_membership; noiseIDs];

end

%% loop: make models and do rejection etc
nBatch = size(Results.Time,2);

ResultsFields = {'normVIQmaxFullOne','normVIConsensusFullOne','normVIQmaxSparseOne','normVIConsensusSparseOne',...
                    'normVIQmaxFullOwn','normVIConsensusFullOwn','normVIQmaxSparseOwn','normVIConsensusSparseOwn',...
                    'normVILouvainOne','normVIMultiwayOne','normVIMultiwayQOne','normVILouvainOwn','normVIMultiwayOwn','normVIMultiwayQOwn','Time',...
                    'nGrpsLouvain','nGrpsMultiway','nGrpsMultiwayQ','nGrpsQmaxFull','nGrpsConsensusFull','nGrpsQmaxSparse','nGrpsConsensusSparse'};
LoopResultsFields = {'QmaxCluster','ConsCluster','normVIQmaxSpectra','normVIConsensusSpectra','nGrpsQmaxSpectra','nGrpsConsensusSpectra'};

fieldsize = [nBatch,1]; % parfor cannot handle matrices of structs...
ClustResults = emptyStruct(ResultsFields,fieldsize);
LoopResults = emptyStruct(LoopResultsFields,fieldsize);

blnP = autoParallel;  % set-up parallel processing with scaled number of cores.

for iF = 1:numel(Model.F_noise)
    % make a temporary variable for broadcast to workers, to prevent
    % machine from seizing up
    % tempNetwork = squeeze(Network(:,iF,:));
    tempNetwork(1:numel(Model.P_of_noise),1:nBatch) = struct('W',[],'ExpW',[],'ixFinal_Sparse',[],'ixFinal_Full',[]);
    % make loop to recast as single
    for iB  = 1:nBatch
        for iP = 1:numel(Model.P_of_noise)
             tempNetwork(iP,iB).W = single(Network(iP,iF,iB).W);
             tempNetwork(iP,iB).ExpW = single(Network(iP,iF,iB).ExpW);
             tempNetwork(iP,iB).ixFinal_Sparse = Network(iP,iF,iB).ixFinal_Sparse;
             tempNetwork(iP,iB).ixFinal_Full = Network(iP,iF,iB).ixFinal_Full;
        end
    end
    
    parfor iB = 1:nBatch
    % for iB = 1:nBatch
    
        for iP = 1:numel(Model.P_of_noise)
            disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_noise)) '; F:' num2str(iF) '/' num2str(numel(Model.F_noise)) '; batch: ' num2str(iB)])

            tic 
            % Louvain on synthetic network
            [LouvCluster,LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(tempNetwork(iP,iB).W,1,1);  % return 1st level of hierarchy only
            ClustResults(iB).nGrpsLouvain(iP,iF) = max(LouvCluster{1}{1});
            [~,ClustResults(iB).normVILouvainOwn(iP,iF)] = VIpartitions(LouvCluster{1}{1},Partition(iF).owngroups);
            [~,ClustResults(iB).normVILouvainOne(iP,iF)] = VIpartitions(LouvCluster{1}{1},Partition(iF).onegroup);
    
            % multi-way spectra on synthetic network
            [bestPartition,maxQPartition] = multiwaySpectCommDet(tempNetwork(iP,iB).W,clusterpars.maxMultiway);
            ClustResults(iB).nGrpsMultiway(iP,iF) = max(bestPartition);
            [~,ClustResults(iB).normVIMultiwayOwn(iP,iF)] = VIpartitions(bestPartition,Partition(iF).owngroups);
            [~,ClustResults(iB).normVIMultiwayOne(iP,iF)] = VIpartitions(bestPartition,Partition(iF).onegroup);

            ClustResults(iB).nGrpsMultiwayQ(iP,iF) = max(maxQPartition);
            [~,ClustResults(iB).normVIMultiwayQOwn(iP,iF)] = VIpartitions(maxQPartition,Partition(iF).owngroups);
            [~,ClustResults(iB).normVIMultiwayQOne(iP,iF)] = VIpartitions(maxQPartition,Partition(iF).onegroup);

            %% clustering on low-D space from rejection: on signal network
            % (1) sparse WCM
           
            Wsignal = tempNetwork(iP,iB).W(tempNetwork(iP,iB).ixFinal_Sparse,tempNetwork(iP,iB).ixFinal_Sparse);        % get final signal matrix
            ExpWsignal = tempNetwork(iP,iB).ExpW(tempNetwork(iP,iB).ixFinal_Sparse,tempNetwork(iP,iB).ixFinal_Sparse);  % extract expected matrix for just final signal
            LoopResults = clusterLowDNetwork(Wsignal,ExpWsignal,Results.SpectraSparseWCM.Groups(iP,iF,iB),Results.SpectraSparseWCM.Groups(iP,iF,iB),clusterpars.nreps,Partition(iF).owngroups(tempNetwork(iP,iB).ixFinal_Sparse));
            [ClustResults(iB).normVIQmaxSparseOwn(iP,iF),ClustResults(iB).normVIConsensusSparseOwn(iP,iF),ClustResults(iB).nGrpsQmaxSparse(iP,iF),ClustResults(iB).nGrpsConsensusSparse(iP,iF)] ...
                    = deal(LoopResults.normVIQmaxSpectra,LoopResults.normVIConsensusSpectra,LoopResults.nGrpsQmaxSpectra,LoopResults.nGrpsConsensusSpectra);
            % also compute one-group VI
            if ~isempty(LoopResults.QmaxCluster)
                ClustResults(iB).normVIQmaxSparseOne(iP,iF) = VIpartitions(LoopResults.QmaxCluster,Partition(iF).onegroup(tempNetwork(iP,iB).ixFinal_Sparse));
            else
                ClustResults(iB).normVIQmaxSparseOne(iP,iF) = 0;
            end
            if ~isempty(LoopResults.ConsCluster)
                ClustResults(iB).normVIConsensusSparseOne(iP,iF) = VIpartitions(LoopResults.ConsCluster,Partition(iF).onegroup(tempNetwork(iP,iB).ixFinal_Sparse));
            else
                ClustResults(iB).normVIConsensusSparseOne(iP,iF) = nan;                
            end
            
            
           
            % (2) full WCM
            Wsignal = tempNetwork(iP,iB).W(tempNetwork(iP,iB).ixFinal_Full,tempNetwork(iP,iB).ixFinal_Full);        % get final signal matrix            
            P = expectedA(tempNetwork(iP,iB).W);
            Psignal = P(tempNetwork(iP,iB).ixFinal_Full,tempNetwork(iP,iB).ixFinal_Full);
            LoopResults = clusterLowDNetwork(Wsignal,Psignal,Results.SpectraFullWCM.Groups(iP,iF,iB),Results.SpectraFullWCM.Groups(iP,iF,iB),clusterpars.nreps,Partition(iF).owngroups(tempNetwork(iP,iB).ixFinal_Full));
            [ClustResults(iB).normVIQmaxFullOwn(iP,iF),ClustResults(iB).normVIConsensusFullOwn(iP,iF),ClustResults(iB).nGrpsQmaxFull(iP,iF),ClustResults(iB).nGrpsConsensusFull(iP,iF)] ...
                    = deal(LoopResults.normVIQmaxSpectra,LoopResults.normVIConsensusSpectra,LoopResults.nGrpsQmaxSpectra,LoopResults.nGrpsConsensusSpectra);
            % also compute one-group VI
            if ~isempty(LoopResults.QmaxCluster)
                ClustResults(iB).normVIQmaxFullOne(iP,iF) = VIpartitions(LoopResults.QmaxCluster,Partition(iF).onegroup(tempNetwork(iP,iB).ixFinal_Full));
                ClustResults(iB).normVIConsensusFullOne(iP,iF) = VIpartitions(LoopResults.ConsCluster,Partition(iF).onegroup(tempNetwork(iP,iB).ixFinal_Full));
            else
                ClustResults(iB).normVIQmaxFullOne(iP,iF) = 0;
                ClustResults(iB).normVIConsensusFullOne(iP,iF) = 0;
                
            end
            
            
            ClustResults(iB).Time(iP,iF) = toc;
        end
    end
end

%% save
save([fpath 'Clustering' fname],'ClustResults','clusterpars')


