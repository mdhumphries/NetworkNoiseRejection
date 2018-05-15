%% script to cluster Full and Signal networks, and compare them
% Mark Humphries, 28/2/2017
clear all; close all;

fname = 'LesMis'; 
blnLabels = 1;      % write node labels? Omit for large networks
fontsize = 6;

clusterpars.nreps = 100;
clusterpars.nLouvain = 5;

% load data
load(['Results/Rejected_' fname])

%% cluster - with noise rejection
% consensus modularity
% [Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,Ngrps,~] = allevsplitConTransitive(Data.Aconnected);

% construct new null model
P = Data.ExpA(Data.ixSignal_Final,Data.ixSignal_Final); % extract relevant part of null model

% or make one
% [Signal.Emodel,~,Vmodel,Signal.ExpA] = RndPoissonConfigModel(Data.Aconnected,pars.N,pars.C,optionsModel);
% P = Data.Aconnected - Signal.ExpA;  % modularity matrix using chosen null model
% % find low-dimensional projection
% [~,~,Signal.Dn,~] = LowDSpace(P,Signal.Emodel,pars.alpha); % to just obtain low-dimensional projection

% then cluster
if Data.Dn > 0
    [Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,ctr] = ...
                                        ConsensusCommunityDetect(Data.Asignal_final,P,1+Data.Dn,1+Data.Dn,clusterpars.nreps);
else
    Connected.QmaxCluster = []; Connected.Qmax = 0; Connected.ConsCluster = []; Connected.ConsQ = 0;
end
% Louvain algorithm
[Connected.LouvCluster,Connected.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.Asignal_final,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only

%% cluster - without noise rejection
if Data.Dn > 0
    [Full.QmaxCluster,Full.Qmax,Full.ConsCluster,Full.ConsQ,~] = ...
                                                ConsensusCommunityDetect(Data.A,Data.ExpA,1+Data.Dn,1+Data.Dn);
else
    Full.QmaxCluster = []; Full.Qmax = 0; Full.ConsCluster = []; Full.ConsQ = 0;
end

[Full.LouvCluster,Full.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.A,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only

%% plot sorted into group order
numConnected = length(Data.ixSignal_Final);

if Data.Dn > 0
    [H,Ix] = plotClusterMap(Data.Asignal_final,Connected.ConsCluster,[],'S');
    title('Consensus clustering')
    plotorder = Data.ixSignal_Final(Ix);
    
    if blnLabels
        % Add node labelss
        set(gca,'Ytick',1:numConnected);
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        % set(gca,'XTickLabelRotation',90);
    end
    % compare to the Qmax solution at the requested number of groups
    [H,Ix] = plotClusterMap(Data.Asignal_final,Connected.QmaxCluster,[],'S');
    title('Qmax clustering')
    if blnLabels
        % Add node labelss
        set(gca,'Ytick',1:numConnected);
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        % set(gca,'XTickLabelRotation',90);
    end
   
end

for i=1:numel(Connected.LouvCluster)
    CLou = Connected.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
    [H,Ix] = plotClusterMap(Data.Asignal_final,CLou,[],'S');
    plotorder = Data.ixSignal_Final(Ix);
    title(['Louvain ' num2str(i)]);
    if blnLabels
        % Add node labels
        set(gca,'Ytick',1:numConnected);
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
    end
    % set(gca,'XTickLabelRotation',90);
    for j = i+1:numel(Connected.LouvCluster)
        CLou2 = Connected.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
        Connected.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log(numConnected);
    end
end

%% without noise rejection
if Data.Dn > 0
    [H,Ix] = plotClusterMap(Data.A,Full.ConsCluster,[],'S');
    title('Consensus clustering of all')
    plotorder = Ix;
    
    if blnLabels
        % Add node labels
        [srt,I] = sort(Full.ConsCluster,'ascend');
        set(gca,'Ytick',1:numel(Data.ixRetain));
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
    end
end

% Louvain algorithm
for i=1:numel(Full.LouvCluster)
    CLou = Full.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
    [HL,Ix] = plotClusterMap(Data.A,CLou,[],'S');
    title(['Full Louvain ' num2str(i)]);
    plotorder = Ix;
    if blnLabels
        % Add node labels
        set(gca,'Ytick',1:numel(Data.ixRetain));
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
    end
    for j = i+1:numel(Full.LouvCluster)
        CLou2 = Full.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
        Full.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log(numel(Data.ixRetain));
    end
end

save(['Results/Clustered_' fname],'Full','Connected','clusterpars')