%% script to assess performance of group detection as tuning network from random to modular 
% Here in script 02: Communities detected by community detection algorithms
% on noise halo
%
%
% Change log:
% 08/06/2018: initial version
% 04/07/2018 : added comparison to ground-truths
% Mark Humphries

clear all; close all;
addpath('../Network_Analysis_Functions/');
addpath('../ZhangNewman2015/');
addpath('../Helper_Functions/')

%% load networks from rejection script
fname = 'P_rejection_SyntheticEqual_Noise_20180611T132723';  % full set of 100 networks per P(within) level
fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';

load([fpath fname])

%% clustering parameters
clusterpars.nreps = 50;
clusterpars.nLouvain = 1;

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

ResultsFields = {'normVIQmaxSpectra','normVIConsensusSpectra','normVILouvain','normVIMultiway','Time',...
                    'nGrpsLouvain','nGrpsMultiway','nGrpsQmaxSpectra','nGrpsConsensusSpectra'};
fieldsize = [nBatch,1]; % parfor cannot handle matrices of structs...
ClustResults = emptyStruct(ResultsFields,fieldsize);

blnP = autoParallel;  % set-up parallel processing with scaled number of cores.

parfor iB = 1:nBatch
% for iB = 1:nBatch
    
    for iP = 1:numel(Model.P_of_noise)
        for iF = 1:numel(Model.F_noise)
            disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_noise)) '; F:' num2str(iF) '/' num2str(numel(Model.F_noise)) '; batch: ' num2str(iB)])

            tic 
            % Louvain on synthetic network
            [LouvCluster,LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Network(iP,iF,iB).W,1,1);  % return 1st level of hierarchy only
            ClustResults(iB).nGrpsLouvain(iP,iF) = max(LouvCluster{1}{1});
            [~,ClustResults(iB).normVILouvainOwn(iP,iF)] = VIpartitions(LouvCluster{1}{1},Partition(iF).owngroups);
            [~,ClustResults(iB).normVILouvainOne(iP,iF)] = VIpartitions(LouvCluster{1}{1},Partition(iF).onegroup);

            % multi-way spectra on synthetic network
            [bestPartition] = multiwaySpectCommDet(Network(iP,iF,iB).W);
            ClustResults(iB).nGrpsMultiway(iP,iF) = max(bestPartition);
            [~,ClustResults(iB).normVIMultiwayOwn(iP,iF)] = VIpartitions(bestPartition,Partition(iF).owngroups);
            [~,ClustResults(iB).normVIMultiwayOne(iP,iF)] = VIpartitions(bestPartition,Partition(iF).onegroup);

             %% specified clustering on synthetic network
            if Results.SpectraWCM.Groups(iP,iF,iB) > 1
                % using final Signal matrix: stripped down to just Signal
                % nodes, and with giant component found and leaves removed
                ExpW = Network(iP,iF,iB).ExpW(Network(iP,iF,iB).ixFinal,Network(iP,iF,iB).ixFinal);  % extract expected matrix for just final signal
                [QmaxCluster,Qmax,ConsCluster,ConsQ,~] = ...
                                                            ConsensusCommunityDetect(Network(iP,iF,iB).WsignalFinal,ExpW,Results.SpectraWCM.Groups(iP,iF,iB),Results.SpectraWCM.Groups(iP,iF,iB),clusterpars.nreps);
                % quality of estimation of retained communities
                if ~isempty(QmaxCluster)
                    ClustResults(iB).nGrpsQmaxSpectra(iP,iF) = max(QmaxCluster);
                    [~,ClustResults(iB).normVIQmaxSpectraOwn(iP,iF)] = VIpartitions(QmaxCluster,Partition(iF).owngroups(Network(iP,iF,iB).ixFinal));  % defined on retained network
                    [~,ClustResults(iB).normVIQmaxSpectraOne(iP,iF)] = VIpartitions(QmaxCluster,Partition(iF).onegroup(Network(iP,iF,iB).ixFinal));  % defined on retained network
                else
                    ClustResults(iB).nGrpsQmaxSpectra(iP,iF) = nan; % error in clustering
                    ClustResults(iB).normVIQmaxSpectraOwn(iP,iF) = 0;  
                    ClustResults(iB).normVIQmaxSpectraOne(iP,iF) = 0;
                    
                end

                if ~isempty(ConsCluster)
                    ClustResults(iB).nGrpsConsensusSpectra(iP,iF) = max(ConsCluster);
                    [~,ClustResults(iB).normVIConsensusSpectraOwn(iP,iF)] = VIpartitions(ConsCluster,Partition(iF).owngroups(Network(iP,iF,iB).ixFinal));  % defined on retained network
                    [~,ClustResults(iB).normVIConsensusSpectraOne(iP,iF)] = VIpartitions(ConsCluster,Partition(iF).onegroup(Network(iP,iF,iB).ixFinal));  % defined on retained network
                    
                else
                    ClustResults(iB).nGrpsConsensusSpectra(iP,iF) = nan; % error in clustering
                    ClustResults(iB).normVIConsensusSpectraOwn(iP,iF) = 0;  
                    ClustResults(iB).normVIConsensusSpectraOne(iP,iF) = 0;      
                end
            else
                % no clusters detected, so set to 0
                QmaxCluster = []; Qmax = 0; ConsCluster = []; ConsQ = 0;
                ClustResults(iB).nGrpsQmaxSpectra(iP,iF) = 0;
                ClustResults(iB).nGrpsConsensusSpectra(iP,iF) = 0;
                ClustResults(iB).normVIConsensusSpectraOwn(iP,iF) = 0;  
                ClustResults(iB).normVIConsensusSpectraOne(iP,iF) = 0;
                ClustResults(iB).normVIQmaxSpectraOwn(iP,iF) = 0;  
                ClustResults(iB).normVIQmaxSpectraOne(iP,iF) = 0;
                
            end
            ClustResults(iB).Time(iP,iF) = toc;
        end
    end
end

%% save
save([fpath 'Clustering' fname],'ClustResults','clusterpars')


