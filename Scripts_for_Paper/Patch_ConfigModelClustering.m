% script to cluster using the answers from the vanilla Weighted Configuration model 
%
% Mark Humphries 09/07/2018

clearvars;
fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';

addpath('../Network_Analysis_Functions/');
addpath('../Helper_Functions/')

% rejection parameters
rejectionpars.N = 100;           % repeats of permutation
rejectionpars.alpha = 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
rejectionpars.C = 1;            % conversion factor for real-valued weights (set=1 for integers)
rejectionpars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

clusterpars.nreps = 50;

blnP = autoParallel;  % set-up parallel processing with scaled number of cores.

%% (1) P(between) = 0.05
fname = 'P_rejection_SyntheticEqual_NoNoise_20180702T124726';
load([fpath fname])  % get W for all networks

nBatch = size(Results.Time,2);
ResultsFields = {'normVIQmaxSpectra','normVIConsensusSpectra','nGrpsQmaxSpectra','nGrpsConsensusSpectra'};
fieldsize = [nBatch,1]; % parfor cannot handle matrices of structs...
ConfigResults = emptyStruct(ResultsFields,fieldsize);




% cluster them



parfor iB = 1:nBatch
% for iB = 1:nBatch
    
    for iP = 1:numel(Model.P_of_within)
        disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_within)) '; batch: ' num2str(iB)])
        
        % rederive WCM results
        ConfigReject = rejectConfigModel(Network(iP,iB).W);
        
        % cluster them 
         if ConfigReject.Dn > 0     % if dimensions exist to cluster in...
            [QmaxCluster,Qmax,ConsCluster,ConsQ,~] = ...
                                                        ConsensusCommunityDetect(Network(iP,iB).W,ConfigReject.P,ConfigReject.Dn+1,ConfigReject.Dn+1,clusterpars.nreps);
            % quality of estimation of retained communities
            if ~isempty(QmaxCluster)
                [~,ConfigResults(iB).normVIQmaxSpectra(iP)]=VIpartitions(QmaxCluster,group_membership);
                ConfigResults(iB).nGrpsQmaxSpectra(iP) = max(QmaxCluster);
            else
                ConfigResults(iB).normVIQmaxSpectra(iP)= 0;
                ConfigResults(iB).nGrpsQmaxSpectra(iP) = nan;                
            end
            
            if ~isempty(ConsCluster)
                [~,ConfigResults(iB).normVIConsensusSpectra(iP)]=VIpartitions(ConsCluster,group_membership);
                ConfigResults(iB).nGrpsConsensusSpectra(iP) = max(ConsCluster);
            else
                ConfigResults(iB).normVIConsensusSpectra(iP)= 0;
                ConfigResults(iB).nGrpsConsensusSpectra(iP) = nan;                
            end
        else
            QmaxCluster = []; Qmax = 0; ConsCluster = []; ConsQ = 0;
            ConfigResults(iB).normVIQmaxSpectra(iP)=0;
            ConfigResults(iB).normVIConsensusSpectra(iP)=0;
            ConfigResults(iB).nGrpsQmaxSpectra(iP) = 0;
            ConfigResults(iB).nGrpsConsensusSpectra(iP) = 0;
        end
    end
end


%% (2) P(between) = 0.15
fname = 'P_rejection_SyntheticEqual_NoNoise_20180705T141726';
load([fpath fname])


%% (3) Unequal groups
fname = 'P_rejection_SyntheticUnequal_NoNoise_20180705T180300';
load([fpath fname])


%% (4) Noise
fname = 'P_rejection_SyntheticEqual_Noise_20180611T132723';
load([fpath fname])
