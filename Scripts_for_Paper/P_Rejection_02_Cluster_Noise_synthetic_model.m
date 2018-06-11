%% script to assess performance of group detection as tuning network from random to modular 
% Here in script 02: Communities detected by community detection algorithms
% on noise halo
%
%
% Change log:
% 8/6/2018: initial version
%
% Mark Humphries

clear all; close all;
addpath('../SyntheticModel/');
addpath('../Network_Spectra_Functions/');
addpath('../ZhangNewman2015/');

%% load networks from rejection script
fname = '...';  % full set of 100 networks per P(within) level
fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';

load([fpath fname])

%% clustering parameters
clusterpars.nreps = 50;
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
    
    for iP = 1:numel(Model.P_of_noise)
        for iF = 1:numel(Model.F_noise)
            disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_noise)) '; F:' num2str(iF) '/' num2str(numel(Model.F_noise)) '; batch: ' num2str(iB)])

            tic 
            % Louvain on synthetic network
            [LouvCluster,LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Network(iP,iF,iB).W,1,1);  % return 1st level of hierarchy only
            ClustResults(iB).nGrpsLouvain(iP,iF) = max(LouvCluster{1}{1});

            % multi-way spectra on synthetic network
            [bestPartition] = multiwaySpectCommDet(Network(iP,iF,iB).W);
            ClustResults(iB).nGrpsMultiway(iP,iF) = max(bestPartition);

             %% specified clustering on synthetic network
            if Results.SpectraWCMGroups(iP,iF,iB) > 1
                % using final Signal matrix: stripped down to just Signal
                % nodes, and with giant component found and leaves removed
                ExpW = Network(iP,iF,iB).ExpW(Network(iP,iF,iB).ixFinal,Network(iP,iF,iB).ixFinal);  % extract expected matrix for just final signal
                [QmaxCluster,Qmax,ConsCluster,ConsQ,~] = ...
                                                            ConsensusCommunityDetect(Network(iP,iF,iB).WsignalFinal,ExpW,Results.SpectraWCMGroups(iP,iF,iB),Results.SpectraWCMGroups(iP,iF,iB),clusterpars.nreps);
                % quality of estimation of retained communities
                if ~isempty(QmaxCluster)
                    ClustResults(iB).nGrpsQmaxSpectra(iP,iF) = max(QmaxCluster);
                else
                    ClustResults(iB).nGrpsQmaxSpectra(iP,iF) = nan;                
                end

                if ~isempty(ConsCluster)
                    ClustResults(iB).nGrpsConsensusSpectra(iP,iF) = max(ConsCluster);
                else
                    ClustResults(iB).nGrpsConsensusSpectra(iP,iF) = nan;                
                end
            else
                QmaxCluster = []; Qmax = 0; ConsCluster = []; ConsQ = 0;
                ClustResults(iB).nGrpsQmaxSpectra(iP,iF) = 0;
                ClustResults(iB).nGrpsConsensusSpectra(iP,iF) = 0;
            end
            ClustResults(iB).Time(iP,iF) = toc;
        end
    end
end

%% save
save([fpath 'Clustering' fname],'ClustResults','clusterpars')


