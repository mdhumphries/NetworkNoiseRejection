%% script to cluster Full and Signal networks, and compare them
%
% 28/2/2017: proper script
% 17/07/2018: added Full WCM clustering
%
% Mark Humphries, 


% uncomment to use as standalone
clearvars; close all;

addpath('../Helper_Functions/')
addpath('../Network_Analysis_Functions/')

% add fname here to use as standalone script
% fname = 'LesMis'; 
fname = 'Allen_Gene_Leaf';
blnPlot = 0;
blnLabels = 1;      % write node labels? Omit for large networks
fontsize = 6;

clusterpars.nreps = 100;
clusterpars.nLouvain = 5;
clusterpars.explore = 'explore';  % allow consensus to use more groups than specified by spectral rejection

% load data
load(['../Results/Rejected_' fname])

%% cluster - with noise rejection


% or make one
% [Signal.Emodel,~,Vmodel,Signal.ExpA] = RndPoissonConfigModel(Data.Aconnected,pars.N,pars.C,optionsModel);
% P = Data.Aconnected - Signal.ExpA;  % modularity matrix using chosen null model
% % find low-dimensional projection
% [~,~,Signal.Dn,~] = LowDSpace(P,Signal.Emodel,pars.alpha); % to just obtain low-dimensional projection

% then cluster: Sparse WCM
if Data.Dn > 0 && numel(Data.ixSignal_Final) > 3
    % construct new null model
    P = Data.ExpA(Data.ixSignal_Final,Data.ixSignal_Final); % extract relevant part of null model
   
    [Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,ctr] = ...
                                        ConsensusCommunityDetect(Data.Asignal_final,P,1+Data.Dn,1+Data.Dn,clusterpars.nreps,[],clusterpars.explore);
    % Louvain algorithm
    [Connected.LouvCluster,Connected.LouvQ,~,~] = LouvainCommunityUDnondeterm(Data.Asignal_final,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only

else
    Connected.QmaxCluster = []; Connected.Qmax = 0; Connected.ConsCluster = []; Connected.ConsQ = 0;
    Connected.LouvCluster = cell(clusterpars.nLouvain,1);
    Connected.LouvQ = cell(clusterpars.nLouvain,1);
end

% full WCM
if Control.Dn > 0 && numel(Control.ixSignal_Final) > 3
    P = Control.P(Control.ixSignal_Final,Control.ixSignal_Final); % extract relevant part of null model

    [Connected.QmaxClusterFullWCM,Connected.QmaxFullWCM,Connected.ConsClusterFullWCM,Connected.ConsQFullWCM,ctr] = ...
                                        ConsensusCommunityDetect(Control.Asignal_final,P,1+Control.Dn,1+Control.Dn,clusterpars.nreps,[],clusterpars.explore);
else
    Connected.QmaxClusterFullWCM = []; Connected.QmaxFullWCM = 0; Connected.ConsClusterFullWCM = []; Connected.ConsQFullWCM = 0;
end


%% cluster - without noise rejection
if Data.Dn > 0
    [Full.QmaxCluster,Full.Qmax,Full.ConsCluster,Full.ConsQ,~] = ...
                                                ConsensusCommunityDetect(Data.A,Data.ExpA,1+Data.Dn,1+Data.Dn,clusterpars.nreps,[],clusterpars.explore);
else
    Full.QmaxCluster = []; Full.Qmax = 0; Full.ConsCluster = []; Full.ConsQ = 0;
end

[Full.LouvCluster,Full.LouvQ,~,~] = LouvainCommunityUDnondeterm(Data.A,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only

%% plot sorted into group order
if blnPlot
    numConnected = length(Data.ixSignal_Final);

    if Data.Dn > 0
        [~,~,Ix] = plotClusterMap(Data.Asignal_final,Connected.ConsCluster,[],[],'S');
        title('Consensus clustering')
        plotorder = Data.ixSignal_Final(Ix);

        if blnLabels
            % Add node labelss
            set(gca,'Ytick',1:numConnected);
            set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
            % set(gca,'XTickLabelRotation',90);
        end
        % compare to the Qmax solution at the requested number of groups
        [~,~,Ix] = plotClusterMap(Data.Asignal_final,Connected.QmaxCluster,[],[],'S');
        title('Qmax clustering')
        if blnLabels
            % Add node labelss
            set(gca,'Ytick',1:numConnected);
            set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
            % set(gca,'XTickLabelRotation',90);
        end


        for i=1:numel(Connected.LouvCluster)
            CLou = Connected.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
            [~,~,Ix] = plotClusterMap(Data.Asignal_final,CLou,[],[],'S');
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
                Connected.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log2(numConnected);
            end
        end
    end
    
    %% without noise rejection

    [~,~,Ix] = plotClusterMap(Data.A,Full.ConsCluster,[],[],'S');
    title('Consensus clustering of all')
    plotorder = Ix;

    if blnLabels
        % Add node labels
        [srt,I] = sort(Full.ConsCluster,'ascend');
        set(gca,'Ytick',1:numel(Data.ixRetain));
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
    end


    % Louvain algorithm
    for i=1:numel(Full.LouvCluster)
        CLou = Full.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
        [~,~,Ix] = plotClusterMap(Data.A,CLou,[],[],'S');
        title(['Full Louvain ' num2str(i)]);
        plotorder = Ix;
        if blnLabels
            % Add node labels
            set(gca,'Ytick',1:numel(Data.ixRetain));
            set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        end
        for j = i+1:numel(Full.LouvCluster)
            CLou2 = Full.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
            Full.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log2(numel(Data.ixRetain));
        end
    end
end

save(['../Results/Clustered_' fname],'Full','Connected','clusterpars')