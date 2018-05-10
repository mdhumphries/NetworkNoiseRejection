%% script to separately visualise clustering results in heatmaps

% Mark Humphries, Mat Evans 28/2/2017

clear all; close all
% fname = 'Allen_Gene_Leaf'; 
% fname = 'LesMis'; 
fname = 'StarWarsNetworkEp5'; 

blnLabels = 1;      % write node labels? Omit for large networks
blnExport = 0;
fontsize = 6;

%% load data
load(['Results/Rejected_' fname],'Data','Rejection')
load(['Results/Clustered_' fname],'Full','Connected')

%% plot clustering of connected signal

numConnected = length(Data.ixSignal_Final);

if Data.Dn > 0
    [H,h,Ix] = plotClusterMap(Data.Asignal_final,Connected.ConsCluster,[],[],'S');
    title('Consensus clustering of signal')
    plotorder = Data.ixSignal_Final(Ix);
    
    if blnLabels
        % Add node labelss
        set(gca,'Ytick',1:numConnected);
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        % set(gca,'XTickLabelRotation',90);
    end
    if blnExport exportPPTfig(gcf,[fname 'Consensus_Signal'],[10 15 7 7]); end

    % compare to the Qmax solution at the requested number of groups
    [H,h,Ix] = plotClusterMap(Data.Asignal_final,Connected.QmaxCluster,[],[],'S');
    plotorder = Data.ixSignal_Final(Ix);
    title('Qmax clustering of signal')
    if blnLabels
        % Add node labelss
        set(gca,'Ytick',1:numConnected);
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        % set(gca,'XTickLabelRotation',90);
    end
   
end

%% plot clustering of full network
if Data.Dn > 0
    [H,h,Ix] = plotClusterMap(Data.A,Full.ConsCluster,[],[],'S');
    title('Consensus clustering of all')
    plotorder = Ix;
    
    if blnLabels
        % Add node labels
        [srt,I] = sort(Full.ConsCluster,'ascend');
        set(gca,'Ytick',1:numel(Data.ixRetain));
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
    end
    
    % compare to the Qmax solution at the requested number of groups
    [H,h,Ix] = plotClusterMap(Data.A,Full.QmaxCluster,[],[],'S');
    plotorder = Ix;
    title('Qmax clustering of all')
    if blnLabels
        % Add node labelss
        set(gca,'Ytick',1:numel(Data.ixRetain));
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        % set(gca,'XTickLabelRotation',90);
    end

end
