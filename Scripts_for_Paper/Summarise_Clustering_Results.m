% gather all results into a Table...
clear all; close all;

fnames = dir('../Results/');

nF = numel(fnames);

netCtr = 0;
for iF = 1:nF
    
    if any(strfind(fnames(iF).name,'Clustered'))
        netCtr = netCtr + 1;
        result(netCtr).NetworkName = fnames(iF).name(11:end-4); % strip out 'Clustered' and .mat
        load(['../Results/' fnames(iF).name]);
        
        % keyboard
        % rejection results on signal network
        result(netCtr).SparseSignal_Consensus_Grps = numel(unique(Connected.ConsCluster)); % number of groups found by consensus clustering exploration
        result(netCtr).FullSignal_Consensus_Grps = numel(unique(Connected.ConsClusterFullWCM)); % number of groups found by consensus clustering exploration
        
        % Louvain results on signal
        if ~isempty(Connected.LouvCluster{1})
            n = cellfun(@(x) numel(unique(x{1})),Connected.LouvCluster);    % number of groups in each Louvain clustering
            result(netCtr).Signal_Louvain_MeanGrps = mean(n); 
            result(netCtr).Signal_Louvain_RangeGrps = range(n); 
            if isfield(Connected,'VI_Louvain')
                ix = find(~tril(ones(size(Connected.VI_Louvain))));
                result(netCtr).Signal_Louvain_MeanVI = mean(Connected.VI_Louvain(ix));
            end
        else
            result(netCtr).Signal_Louvain_MeanGrps = 1;
            result(netCtr).Signal_Louvain_RangeGrps = 0;
            result(netCtr).Signal_Louvain_MeanVI = 0;
        end
        
        % results on raw network
        result(netCtr).Raw_Consensus_Grps = numel(unique(Full.ConsCluster)); 
        n = cellfun(@(x) numel(unique(x{1})),Full.LouvCluster);    % number of groups in each Louvain clustering
        result(netCtr).Raw_Louvain_MeanGrps = mean(n); 
        result(netCtr).Raw_Louvain_RangeGrps = range(n); 
        if isfield(Full,'VI_Louvain')
            ix = find(~tril(ones(size(Full.VI_Louvain))));
            result(netCtr).Raw_Louvain_MeanVI = mean(Full.VI_Louvain(ix));
        end
    end
end

Network_Clustering_Table = struct2table(result);
save('../Results/Network_Clustering_Table','Network_Clustering_Table');