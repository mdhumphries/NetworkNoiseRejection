% script to cluster using the answers from the vanilla Weighted Configuration model 
%
% (1) Rederive the spectral rejection from sampling config model
% (1b) do node rejection if doing noise networks...
% (2) Cluster it
%   
% Mark Humphries 09/07/2018

clearvars;
fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';

addpath('../Network_Analysis_Functions/');
addpath('../Network_Spectra_Functions/');
addpath('../Helper_Functions/')

% rejection parameters
rejectionpars.N = 100;           % repeats of permutation
rejectionpars.alpha = 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
rejectionpars.C = 1;            % conversion factor for real-valued weights (set=1 for integers)
rejectionpars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

clusterpars.nreps = 50;

% NodeRejection options
optionsReject.Weight = 'linear'; % 'linear' is default
optionsReject.Norm = 'L2';       % L2 is default

blnP = autoParallel;  % set-up parallel processing with scaled number of cores.
ResultsFields = {'normVIQmaxSpectra','normVIConsensusSpectra','nGrpsQmaxSpectra','nGrpsConsensusSpectra'};

%% (1) P(between) = 0.05
fname = 'P_rejection_SyntheticEqual_NoNoise_20180702T124726';
load([fpath fname])  % get W for all networks

nBatch = size(Results.Time,2);
fieldsize = [nBatch,1]; % parfor cannot handle matrices of structs...
ConfigResults = emptyStruct(ResultsFields,fieldsize);
LoopResults = emptyStruct(ResultsFields,fieldsize);

% parfor iB = 1:nBatch
for iB = 1:nBatch
    
    for iP = 1:numel(Model.P_of_within)
        disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_within)) '; batch: ' num2str(iB)])
        
        % rederive WCM results
        ConfigReject = rejectConfigModel(Network(iP,iB).W,rejectionpars);  % uses parallel toolbox internally
        
        nGrps = ConfigReject.Dn+1;
        % cluster them
        LoopResults = clusterLowDNetwork(Network(iP,iB).W,ConfigReject.P,nGrps,nGrps,clusterpars.nreps,group_membership);
        [ConfigResults(iB).normVIQmaxSpectra(iP),ConfigResults(iB).normVIConsensusSpectra(iP),ConfigResults(iB).nGrpsQmaxSpectra(iP),ConfigResults(iB).nGrpsConsensusSpectra(iP)] ...
                = deal(LoopResults.normVIQmaxSpectra,LoopResults.normVIConsensusSpectra,LoopResults.nGrpsQmaxSpectra,LoopResults.nGrpsConsensusSpectra);
    end
end
        
save([fpath fname '_ClusterConfig'],'ConfigResults');

%% (2) P(between) = 0.15
fname = 'P_rejection_SyntheticEqual_NoNoise_20180705T141726';
load([fpath fname])

nBatch = size(Results.Time,2);
fieldsize = [nBatch,1]; % parfor cannot handle matrices of structs...
ConfigResults = emptyStruct(ResultsFields,fieldsize);
LoopResults = emptyStruct(ResultsFields,fieldsize);

% parfor iB = 1:nBatch
for iB = 1:nBatch
    
    for iP = 1:numel(Model.P_of_within)
        disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_within)) '; batch: ' num2str(iB)])
        
        % rederive WCM results
        ConfigReject = rejectConfigModel(Network(iP,iB).W,rejectionpars);  % uses parallel toolbox internally
        
        nGrps = ConfigReject.Dn+1;
        % cluster them
        LoopResults = clusterLowDNetwork(Network(iP,iB).W,ConfigReject.P,nGrps,nGrps,clusterpars.nreps,group_membership);
        [ConfigResults(iB).normVIQmaxSpectra(iP),ConfigResults(iB).normVIConsensusSpectra(iP),ConfigResults(iB).nGrpsQmaxSpectra(iP),ConfigResults(iB).nGrpsConsensusSpectra(iP)] ...
                = deal(LoopResults.normVIQmaxSpectra,LoopResults.normVIConsensusSpectra,LoopResults.nGrpsQmaxSpectra,LoopResults.nGrpsConsensusSpectra);
    end
end
        
save([fpath fname '_ClusterConfig'],'ConfigResults');


%% (3) Unequal groups
fname = 'P_rejection_SyntheticUnequal_NoNoise_20180705T180300';
load([fpath fname])

nBatch = size(Results.Time,2);
fieldsize = [nBatch,1]; % parfor cannot handle matrices of structs...
ConfigResults = emptyStruct(ResultsFields,fieldsize);
LoopResults = emptyStruct(ResultsFields,fieldsize);

% parfor iB = 1:nBatch
for iB = 1:nBatch
    
    for iP = 1:numel(Model.P_of_within)
        disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_within)) '; batch: ' num2str(iB)])
        
        % rederive WCM results
        ConfigReject = rejectConfigModel(Network(iP,iB).W,rejectionpars);  % uses parallel toolbox internally
        
        nGrps = ConfigReject.Dn+1;
        % cluster them
        LoopResults = clusterLowDNetwork(Network(iP,iB).W,ConfigReject.P,nGrps,nGrps,clusterpars.nreps,group_membership);
        [ConfigResults(iB).normVIQmaxSpectra(iP),ConfigResults(iB).normVIConsensusSpectra(iP),ConfigResults(iB).nGrpsQmaxSpectra(iP),ConfigResults(iB).nGrpsConsensusSpectra(iP)] ...
                = deal(LoopResults.normVIQmaxSpectra,LoopResults.normVIConsensusSpectra,LoopResults.nGrpsQmaxSpectra,LoopResults.nGrpsConsensusSpectra);
    end
end
        
save([fpath fname '_ClusterConfig'],'ConfigResults');


%% (4) Noise
fname = 'P_rejection_SyntheticEqual_Noise_20180611T132723';
load([fpath fname])
NoiseResultsFields = {'normVIQmaxSpectraOwn','normVIConsensusSpectraOwn','normVIQmaxSpectraOne','normVIConsensusSpectraOne','nGrpsQmaxSpectra','nGrpsConsensusSpectra'};
LoopResultsFields = {'QmaxCluster','ConsCluster','normVIQmaxSpectra','normVIConsensusSpectra','nGrpsQmaxSpectra','nGrpsConsensusSpectra'};

nBatch = size(Results.Time,2);
fieldsize = [nBatch,1]; % parfor cannot handle matrices of structs...

ConfigResults = emptyStruct(NoiseResultsFields,fieldsize);
LoopResults = emptyStruct(LoopResultsFields,fieldsize);

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

for iF = 1:numel(Model.F_noise)
    % make a temporary variable for broadcast to workers, to prevent
    % machine from seizing up if we do parfor over batch
    tempNetwork = squeeze(Network(:,iF,:));

    % parfor iB = 1:nBatch
    for iB = 1:nBatch

        for iP = 1:numel(Model.P_of_noise)
            disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_noise)) '; F:' num2str(iF) '/' num2str(numel(Model.F_noise)) '; batch: ' num2str(iB)])

            % rederive WCM results, now with node rejection too
            ConfigReject = rejectConfigModel(tempNetwork(iP,iB).W,rejectionpars,'Reject',optionsReject);  % uses parallel toolbox internally
            
            nGrps = ConfigReject.Dn+1;
                        
            % cluster the signal matrix
            Psignal = ConfigReject.P(ConfigReject.ixSignal_Final,ConfigReject.ixSignal_Final);
            
            LoopResults = clusterLowDNetwork(ConfigReject.Asignal_final,Psignal,nGrps,nGrps,clusterpars.nreps,Partition(iF).owngroups(ConfigReject.ixSignal_Final));
            [ConfigResults(iB).normVIQmaxSpectraOwn(iP,iF),ConfigResults(iB).normVIConsensusSpectraOwn(iP,iF),ConfigResults(iB).nGrpsQmaxSpectra(iP,iF),ConfigResults(iB).nGrpsConsensusSpectra(iP,iF)] ...
                    = deal(LoopResults.normVIQmaxSpectra,LoopResults.normVIConsensusSpectra,LoopResults.nGrpsQmaxSpectra,LoopResults.nGrpsConsensusSpectra);
            % now check against one-group ground truth
            if ~isempty(LoopResults.QmaxCluster)
                ConfigResults(iB).normVIQmaxSpectraOne(iP,iF) = VIpartitions(LoopResults.QmaxCluster,Partition(iF).onegroup(ConfigReject.ixSignal_Final));
                ConfigResults(iB).normVIConsensusSpectraOne(iP,iF) = VIpartitions(LoopResults.ConsCluster,Partition(iF).onegroup(ConfigReject.ixSignal_Final));
            else
                ConfigResults(iB).normVIQmaxSpectraOne(iP,iF) = 0;
                ConfigResults(iB).normVIConsensusSpectraOne(iP,iF) = 0;
                
            end
        end
    end
end

save([fpath fname '_ClusterConfig'],'ConfigResults','Partition');

