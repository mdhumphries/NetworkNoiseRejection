% Script to run analysis of planted noise
% Abhinav Singh & Mark Humphries

clear all; close all

addpath('ZhangNewman2015/');

%% params for WSBM model
%static
cc=[100, 20];       % within communities


%looping
fraction_periphery_grid=[0.2 0.5 0.8 1.0];
gg_cc_grid=[0.1 0.3 0.5 0.7 0.9];   % between communities: proportion of within communities' weight parameters
cp_gg_grid=[0.1 0.3 0.5 0.7 0.9 1];   % from periphery to communities: proportion of between communities' weight parameters

%% params for noise rejection

pars.N = 100;           % repeats of permutation
pars.alpha = 0; %0.95; % 0.95; % 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
pars.Model = 'Poiss';   % or 'WCM' . % which null model
pars.C = 1;             % conversion factor for real-valued weights (set=1 for integers)
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
for count_peri=1:numel(fraction_periphery_grid)
    fraction_periphery=fraction_periphery_grid(count_peri);
    for count_gg_cc=1:numel(gg_cc_grid)
        gg_cc=gg_cc_grid(count_gg_cc);
        mixing_index=gg_cc;
        for count_cp_gg=1:numel(cp_gg_grid)
            count_cp_gg
            cp_gg=cp_gg_grid(count_cp_gg);
            noise_index=cp_gg;
            mixing=mixing_index * cc;
            noise=noise_index * mixing; 
            pp=noise/2;
            %generate planted noise network
            options.weight_dist=struct('ingroup', cc, 'outgroup', mixing,...
                'periphery_periphery', pp, 'periphery_core', noise);
            network=test_noise_rejection_planted_noise(50,3, 'low',fraction_periphery, options);
            % unit analyses
            
            %%% HERE: replace with function based on main template
            %%% script....
            Data = reject_the_noise(network.adjacency,pars,optionsModel,optionsReject);
            
            true_members=find(network.membership >= 0);
            performance=numel(intersect(true_members,Data.ixSignal_Final))/numel(true_members);  % proportion of true members retained
            
            communites_estimated_core=network.membership(Data.ixSignal_Final) +1;  % the community membership of the retained nodes
            
                       
            %% do all 3 detection methods
            % construct new null mode
            P = Data.ExpA(Data.ixSignal_Final,Data.ixSignal_Final); % extract relevant part of null model

            % then cluster
            if Data.Dn > 0
                [Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,ctr] = ...
                                                    ConsensusCommunityDetect(Data.Asignal_final,P,1+Data.Dn,1+Data.Dn,clusterpars.nreps);
            else
                Connected.QmaxCluster = []; Connected.Qmax = 0; Connected.ConsCluster = []; Connected.ConsQ = 0;
            end
            % quality of estimation of retained communities
            nVIQmaxSpectra=VIpartitions(Connected.QmaxCluster,communites_estimated_core) / log(numel(network.membership));
            nVIConsensusSpectra=VIpartitions(Connected.ConsCluster,communites_estimated_core) / log(numel(network.membership));
            
            % Louvain algorithm
            [Connected.LouvCluster,Connected.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.Asignal_final,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
            for j = 1:clusterpars.nLouvain
                CLou = Connected.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
                Connected.VI_Louvain(j) = VIpartitions(CLou,communites_estimated_core) ./ log(numel(network.membership));
            end
            nVILouvainMin = min(Connected.VI_Louvain);
            nVILouvainMax = max(Connected.VI_Louvain);
            
            % multi-way spectra
            [bestPartition] = multiwaySpectCommDet(Data.Asignal_final);
            nVIMultiway = VIpartitions(bestPartition,communites_estimated_core) / log(numel(network.membership));
                     
            temp_table=table(fraction_periphery, mixing_index, noise_index, performance, nVIQmaxSpectra, nVIConsensusSpectra,nVILouvainMin,nVILouvainMax,nVIMultiway);
            core_perf=[core_perf; temp_table];
            
            
            
        end
    end

end
%% visualise loop analyses

fname = datestr(now,30);

save(['Results/WSBM_' fname],'core_perf','pars','optionsModel','optionsReject','clusterpars')

