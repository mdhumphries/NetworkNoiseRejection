%% script to cluster Full and Signal networks, and compare them
% Mark Humphries, 28/2/2017
clear all; close all;

fname = 'Lesmis';
nLouvain = 5;
fontsize = 6;

%% load data
load(['Results/Rejected_' fname],'Data')

%% cluster - with noise rejection
% consensus modularity
% [C,Qmax,Ccon,Qc,Ngrps,Q] = allevsplitConTransitive(Asignal);
[Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,Ngrps,~] = allevsplitConTransitive(Data.Aconnected);

% Louvain algorithm
[Connected.LouvCluster,Connected.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.Aconnected,nLouvain,1);  % run 5 times; return 1st level of hierarchy only

%% cluster - without noise rejection
[Full.QmaxCluster,Full.Qmax,Full.ConsCluster,Full.ConsQ,Ngrps,~] = allevsplitConTransitive(Data.A);

[Full.LouvCluster,Full.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.A,nLouvain,1);  % run 5 times; return 1st level of hierarchy only

%% plot sorted into group order

[H,Ix] = plotClusterMap(Data.Aconnected,Connected.ConsCluster,[],'S');
title('Consensus clustering')
plotorder = Data.ixConnectedSignal(Ix);

% Add node labels
numConnected = length(Data.ixConnectedSignal);
set(gca,'Ytick',1:numConnected);
set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
% set(gca,'XTickLabelRotation',90);

for i=1:numel(Connected.LouvCluster)
    CLou = Connected.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
    [H,Ix] = plotClusterMap(Data.Aconnected,CLou,[],'S');
    plotorder = Data.ixConnectedSignal(Ix);
    title(['Louvain ' num2str(i)]);
    % Add node labels
    set(gca,'Ytick',1:numConnected);
    set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
    % set(gca,'XTickLabelRotation',90);
    for j = i+1:numel(Connected.LouvCluster)
        CLou2 = Connected.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
        Connected.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log(numConnected);
    end
end

%% without noise rejection

[H,Ix] = plotClusterMap(Data.A,Full.ConsCluster,[],'S');
title('Consensus clustering of all')
plotorder = Ix;

% Add node labels
[srt,I] = sort(Full.ConsCluster,'ascend');
set(gca,'Ytick',1:numel(Data.ixRetain));
set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);


% Louvain algorithm
for i=1:numel(Full.LouvCluster)
    CLou = Full.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
    [HL,Ix] = plotClusterMap(Data.A,CLou,[],'S');
    title(['Full Louvain ' num2str(i)]);
    plotorder = Ix;
    % Add node labels
    % [srt,I] = sort(CLou,'ascend');
    set(gca,'Ytick',1:numel(Data.ixRetain));
    set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
    for j = i+1:numel(Full.LouvCluster)
        CLou2 = Full.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
        Full.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log(numel(Data.ixRetain));
    end
end

save(['Results/Clustered_' fname],'Full','Connected')