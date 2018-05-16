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
        result(netCtr).Signal_Consensus_Grps = numel(unique(Connected.ConsCluster)); 
        n = cellfun(@(x) numel(unique(x{1})),Connected.LouvCluster);    % number of groups in each Louvain clustering
        result(netCtr).Signal_Louvain_MeanGrps = mean(n); 
        result(netCtr).Signal_Louvain_RangeGrps = range(n); 
        if isfield(Connected,'VI_Louvain')
            ix = find(~tril(ones(size(Connected.VI_Louvain))));
            result(netCtr).Signal_Louvain_MeanVI = mean(Connected.VI_Louvain(ix));
        end
        result(netCtr).Full_Consensus_Grps = numel(unique(Full.ConsCluster)); 
        n = cellfun(@(x) numel(unique(x{1})),Full.LouvCluster);    % number of groups in each Louvain clustering
        result(netCtr).Full_Louvain_MeanGrps = mean(n); 
        result(netCtr).Full_Louvain_RangeGrps = range(n); 
        if isfield(Full,'VI_Louvain')
            ix = find(~tril(ones(size(Full.VI_Louvain))));
            result(netCtr).Full_Louvain_MeanVI = mean(Full.VI_Louvain(ix));
        end
    end
end

Network_Clustering_Table = struct2table(result);
save('../Results/Network_Clustering_Table','Network_Clustering_Table');