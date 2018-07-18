%% script to reject and cluster all selected data networks
% 
% Mark Humphries 17/07/2018

clearvars;

networks = {'StarWarsNetworkEp1','StarWarsNetworkEp2','StarWarsNetworkEp3','StarWarsNetworkEp4','StarWarsNetworkEp5','StarWarsNetworkEp6',...
            'Allen_Gene_Leaf','Lesmis','adjnoun','polblogs','dolphins','cosyneFinalData','celegansneural','power'};

nNet = numel(networks);
% sanity check: run through list, and check all can be loaded without error
for iN = 1:nNet
    load(['../Networks/' networks{iN}]);
end

%% analyse each network
for iN = 1:nNet
    fname = networks{iN}
    
%     % run rejection
%     clearvars -except fname networks iN nNet
%     run Data_Noise_Rejection_Script.m
    
    clearvars -except fname networks iN nNet
    % run clustering
    run Cluster_Network.m
end