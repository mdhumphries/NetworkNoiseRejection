%% script to assess performance on equal-sized groups 
% 
% Mark Humphries 7/8/2017

clear all; close all;
addpath('SyntheticModel/');
addpath('ZhangNewman2015/');

%% fixed parameters edge parameters
Model.N = [100,100,100];  % size of modules

Model.P.between = 0.05;   

% stength distribution parameters
Model.Spar.distribution =  'Poisson';  % type of distribution
Model.Spar.a = 50;                    % scale: in addition to existing edges
Model.Spar.b = 1;                     % spread

%% range of parameters
Model.P_of_within = [0.05:0.025:0.2];
Model.alpha_range = [0.5:0.025:0.6 0.65:0.05:0.8];

%% rejection and clustering parameters
rejectionpars.N = 100;           % repeats of permutation
rejectionpars.alpha = 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
rejectionpars.Model = 'Poiss';   % or 'WCM' . % which null model
rejectionpars.C = 1;            % conversion factor for real-valued weights (set=1 for integers)
rejectionpars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble? 
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?

% NodeRejection options
optionsReject.Weight = 'linear'; % 'linear' is default
optionsReject.Norm = 'L2';       % L2 is default

clusterpars.nreps = 100;
clusterpars.nLouvain = 5;

options.blnViz = 0;


%% loop: make models and do rejection etc
n = sum(Model.N);
P.between = Model.P.between;

% group membership
Nsum = [0 cumsum(Model.N)];
group_membership = zeros(n,1);
for iG = 1:numel(Nsum)-1
    group_membership(Nsum(iG)+1:Nsum(iG+1)) = iG;
end 

% later: loop here to resample degree sequence
Model.S = sample_strength(n,Model.Spar); % sample strength distribution according to current rules

for iP = 1:numel(Model.P_of_within)
    for iA = 1:numel(Model.alpha_range)
        iA
        tic
        % assign parameters
        
        P.in = Model.P_of_within(iP);
        alpha = Model.alpha_range(iA);
        
        %% make model
        A = wire_edges(Model.N,P);  % make adjacency matrix
        W = weight_edges(A,Model.N,Model.S,alpha); % use Poisson generative model to create Weight matrix
        
        %% do spectral rejection on model
        [Data,Rejection,Control] = reject_the_noise(W,rejectionpars,optionsModel,optionsReject);        

        % RESULTS:
        % number of groups recovered by spectra vs same count from other
        % approaches
        Results.SpectraWCMGroups(iP,iA) = Data.Dn+1;
        Results.SpectraConfigGroups(iP,iA) = Control.Dn+1;
        Results.PosEigWCMGroups(iP,iA) = Data.PosDn+1;
        Results.PosEigConfigGroups(iP,iA) = Control.PosDn+1;
        
        % keyboard           
        %% do all 3 detection methods on Signal
        % construct new null model

        %% and all again, on full matrix, assigning noise to its own community
        if Data.Dn > 0
            [Full.QmaxCluster,Full.Qmax,Full.ConsCluster,Full.ConsQ,~] = ...
                                                        ConsensusCommunityDetect(Data.A,Data.ExpA,1+Data.Dn,1+Data.Dn);
        else
            Full.QmaxCluster = []; Full.Qmax = 0; Full.ConsCluster = []; Full.ConsQ = 0;
        end
        % quality of estimation of retained communities
        Results.nVIFull_QmaxSpectra(iP,iA)=VIpartitions(Full.QmaxCluster,group_membership) / log(numel(group_membership));
        Results.nVIFull_ConsensusSpectra(iP,iA)=VIpartitions(Full.ConsCluster,group_membership) / log(numel(group_membership));

        [Full.LouvCluster,Full.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.A,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
        for j = 1:clusterpars.nLouvain
            CLou = Full.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
            Full.VI_Louvain(j) = VIpartitions(CLou,group_membership) ./ log(numel(group_membership));
        end
        Results.nVIFull_LouvainMin(iP,iA) = min(Full.VI_Louvain);
        Results.nVIFull_LouvainMax(iP,iA) = max(Full.VI_Louvain);

        % multi-way spectra
        [bestPartition] = multiwaySpectCommDet(Data.A);
        Results.nVIFull_Multiway(iP,iA) = VIpartitions(bestPartition,group_membership) / log(numel(group_membership));
        
        Results.Time(iP,iA) = toc
    end
end

fname = datestr(now,30);

save(['Results/SyntheticEqual_NoNoise_' fname],'Results','Model','rejectionpars','optionsModel','optionsReject','clusterpars')


