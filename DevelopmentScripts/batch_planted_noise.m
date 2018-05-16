% Script to run analysis of planted noise
% Abhinav Singh & Mark Humphries

clear all; close all

addpath('ZhangNewman2015/');

%% params for WSBM model
%% communities
Communities.size = 50;
Communities.N = 3;

%static
Communities.cc=[100, 20];       % within communities


%looping
Communities.fraction_periphery_grid=[0.2 0.5 0.8 1.0];
Communities.gg_cc_grid=[0.1 0.3 0.5 0.7 0.9];   % between communities: proportion of within communities' weight parameters
Communities.cp_gg_grid=[0.1 0.3 0.5 0.7 0.9 1];   % from periphery to communities: proportion of between communities' weight parameters

T  =numel(Communities.fraction_periphery_grid) * numel(Communities.gg_cc_grid) * numel(Communities.cp_gg_grid);
%% params for noise rejection

pars.N = 100;           % repeats of permutation
pars.alpha = 0.95; %0.95; % 0.95; % 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
pars.Model = 'Poiss';   % or 'WCM' . % which null model
pars.C = 100;             % conversion factor for real-valued weights (set=1 for integers)
pars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble? 
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?

% NodeRejection options
optionsReject.Weight = 'linear'; % 'linear' is default
optionsReject.Norm = 'L2';       % L2 is default

clusterpars.nreps = 100;
clusterpars.nLouvain = 5;

options.blnViz = 0;

%% variables
core_perf=table;


%% set loops for params
for count_peri=1:numel(Communities.fraction_periphery_grid)
    count_peri
    fraction_periphery=Communities.fraction_periphery_grid(count_peri);
    for count_gg_cc=1:numel(Communities.gg_cc_grid)
        gg_cc=Communities.gg_cc_grid(count_gg_cc);
        mixing_index=gg_cc;
        for count_cp_gg=1:numel(Communities.cp_gg_grid)
            count_cp_gg
            tic
            cp_gg=Communities.cp_gg_grid(count_cp_gg);
            noise_index=cp_gg;
            mixing=mixing_index * Communities.cc;
            noise=noise_index * mixing; 
            pp=noise/2;
            %generate planted noise network
            options.weight_dist=struct('ingroup', Communities.cc, 'outgroup', mixing,...
                'periphery_periphery', pp, 'periphery_core', noise);
            network=test_noise_rejection_planted_noise(Communities.size,Communities.N, 'low',fraction_periphery, options);
            % unit analyses
            
            %%% Main dimension reduction and node rejection function
            [Data,Rejection] = reject_the_noise(network.adjacency,pars,optionsModel,optionsReject);        
            
            true_members=find(network.membership >= 0);
            Nperiph = numel(network.membership) - numel(true_members);
            signal=numel(intersect(true_members,Data.ixSignal_Final))/numel(true_members);  % proportion of true members retained
            additional = numel(setdiff(Data.ixSignal_Final,true_members))/Nperiph; % proportion of periphery retained
            
            % keyboard           
            %% do all 3 detection methods on Signal
            % construct new null mode
            P = Data.ExpA(Data.ixSignal_Final,Data.ixSignal_Final); % extract relevant part of null model
            
            % set of comparison communities here
            communites_estimated_core=network.membership(Data.ixSignal_Final) +1;  % the community membership of the retained nodes
            
            % then cluster
            if Data.Dn > 0
                [Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,ctr] = ...
                                                    ConsensusCommunityDetect(Data.Asignal_final,P,1+Data.Dn,1+Data.Dn,clusterpars.nreps);
            else
                Connected.QmaxCluster = []; Connected.Qmax = 0; Connected.ConsCluster = []; Connected.ConsQ = 0;
            end
            
            % quality of estimation of retained communities
            nVIQmaxSpectra=VIpartitions(Connected.QmaxCluster,communites_estimated_core) / log(numel(communites_estimated_core));
            nVIConsensusSpectra=VIpartitions(Connected.ConsCluster,communites_estimated_core) / log(numel(communites_estimated_core));
            
            % Louvain algorithm
            [Connected.LouvCluster,Connected.LouvQ,~,~] = LouvainCommunityUDnondeterm(Data.Asignal_final,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
            for j = 1:clusterpars.nLouvain
                CLou = Connected.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
                Connected.VI_Louvain(j) = VIpartitions(CLou,communites_estimated_core) ./ log(numel(communites_estimated_core));
            end
            nVILouvainMin = min(Connected.VI_Louvain);
            nVILouvainMax = max(Connected.VI_Louvain);
            
            % multi-way spectra
            [bestPartition] = multiwaySpectCommDet(Data.Asignal_final);
            nVIMultiway = VIpartitions(bestPartition,communites_estimated_core) / log(numel(communites_estimated_core));
            
            %% and all again, on full matrix, assigning noise to its own community
            if Data.Dn > 0
                [Full.QmaxCluster,Full.Qmax,Full.ConsCluster,Full.ConsQ,~] = ...
                                                            ConsensusCommunityDetect(Data.A,Data.ExpA,1+Data.Dn,1+Data.Dn);
            else
                Full.QmaxCluster = []; Full.Qmax = 0; Full.ConsCluster = []; Full.ConsQ = 0;
            end
            % quality of estimation of retained communities
            nVIFull_QmaxSpectra=VIpartitions(Full.QmaxCluster,network.membership+1) / log(numel(network.membership));
            nVIFull_ConsensusSpectra=VIpartitions(Full.ConsCluster,network.membership+1) / log(numel(network.membership));

            [Full.LouvCluster,Full.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.A,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
            for j = 1:clusterpars.nLouvain
                CLou = Full.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
                Full.VI_Louvain(j) = VIpartitions(CLou,network.membership+1) ./ log(numel(network.membership));
            end
            nVIFull_LouvainMin = min(Full.VI_Louvain);
            nVIFull_LouvainMax = max(Full.VI_Louvain);

            % multi-way spectra
            [bestPartition] = multiwaySpectCommDet(Data.A);
            nVIFull_Multiway = VIpartitions(bestPartition,network.membership+1) / log(numel(network.membership));
            
            %% all results in a table
            temp_table=table(fraction_periphery, mixing_index, noise_index, signal, additional, ...
                                nVIQmaxSpectra, nVIConsensusSpectra,nVILouvainMin,nVILouvainMax,nVIMultiway,...
                                nVIFull_QmaxSpectra, nVIFull_ConsensusSpectra,nVIFull_LouvainMin,nVIFull_LouvainMax,nVIFull_Multiway);
            core_perf=[core_perf; temp_table];
            
            
            toc
        end
    end

end
%% visualise loop analyses

fname = datestr(now,30);

save(['Results/WSBM_' fname],'core_perf','Communities','pars','optionsModel','optionsReject','clusterpars')

