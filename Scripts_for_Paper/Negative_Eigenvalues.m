%% negative eigenvalues
%
% Mark Humphries 26/7/2018

clearvars; close all

load('../Results/Network_Rejection_Table')

% location of null model eigenvectors
if ispc
    storepath = 'C:\Users\lpzmdh\Dropbox\Analyses\Networks\DataNets_Null_EigVectors\';
else
    storepath = '/Users/mqbssmhg/Dropbox/Analyses/Networks/DataNets_Null_EigVectors/';
end

addpath('../Network_Spectra_Functions/')
addpath('../Network_Analysis_Functions/')
addpath('../Helper_Functions/')

fontsize = 6;
nreps = 50;     % of each distance metric
dims = 'all';   % use all embedding dimensions for each k-means clustering

%% find networks with negative eigenvalues below lower limit
fnames = Network_Rejection_Table.NetworkName(Network_Rejection_Table.SparseWCM_NegDn > 0);
NNegDims = Network_Rejection_Table.SparseWCM_NegDn(Network_Rejection_Table.SparseWCM_NegDn > 0);

for iF = 1:numel(fnames)
    % get rejection results
    load(['../Results/Rejected_' fnames{iF}])
    
    % plot eigenvector
    [Nvec,ix] = sort(Data.Nspace(:,1));
    
    % ytick = find(abs(Nvec) > prctile(abs(Nvec),90));
    ytick = [1 numel(Nvec)];
    figure
    barh(Nvec)
    ylabel('Weight')
    set(gca,'YTick',ytick,'YTickLabel',Data.nodelabels(ix(ytick),:))
    title(fnames{iF})
    exportPPTfig(gcf,[fnames{iF} '_EigVec1'],[10 15 9 9])
    
    if NNegDims(iF) > 1
        [Nvec,ix] = sort(Data.Nspace(:,2));
        % ytick = find(abs(Nvec) > prctile(abs(Nvec),90));
        tick = [1 numel(Nvec)];
        figure
        barh(Nvec)
        ylabel('Weight')
        set(gca,'YTick',ytick,'YTickLabel',Data.nodelabels(ix(ytick),:))
        title([fnames{iF} ' 2nd negative eigenvector'])   
        exportPPTfig(gcf,[fnames{iF} '_EigVec2'],[10 15 9 9])
    end
    
    %% load null model eigenvectors
    load([storepath '/NullModel_Eigenspectrum_' fnames{iF}])
    
    %% do node rejection in this negative-D space
    B = Data.A - Data.ExpA;  % modularity matrix using chosen null model
    % NB parameters have been loaded alongside Data
    
    % use rejection parts, but on lower-bound...
    optionsReject.Bounds = 'Lower';
    Rejection = NodeRejection(B,Data.Emodel,pars.I,Vmodel,optionsReject); % N.B. also calls LowDSpace function to find projections   
    
    % proportions rejected
    PropNodesRetained(iF) = numel(Rejection.ixSignal) ./ numel(Nvec);
    
    %% make into clusters
    Asignal = Data.A(Rejection.ixSignal,Rejection.ixSignal); 

    if NNegDims(iF) == 1
        % if 1 eigenvector, just +/- entries after node rejection
        Grp_Neg = ones(numel(Rejection.ixSignal),1);
        Grp_Neg(Data.Nspace(Rejection.ixSignal) > 0) = 2;

        [H,C,I] = plotClusterMap(Asignal,Grp_Neg,[],[],'S'); 
         plotorder = Rejection.ixSignal(I);
    
        % Add node labels
        set(gca,'Ytick',1:length(Rejection.ixSignal));
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        % set(gca,'XTickLabelRotation',90);
        title(fnames{iF})
        exportPPTfig(gcf,[fnames{iF} '_ClusterMap'],[10 15 12 12])
    else
        % keyboard
    %% if 2 or more, pass to Qmax....
    
        % refine signal matrix for clustering
        % connected signal matrix: find largest component, and use that - store
        % others
        [Asignal_comp,ixRetain,~,~] = prep_A(Asignal); 
        ixSignal_comp = Rejection.ixSignal(ixRetain);  % original node indices

        % and then strip out leaves - nodes with single links
        K = sum(Asignal_comp);
        ixLeaves = find(K==1); ixKeep = find(K > 1);

        ixSignal_Final = ixSignal_comp(ixKeep);
        ixSignal_Leaves = ixSignal_comp(ixLeaves);
        Asignal_final = Asignal_comp(ixKeep,ixKeep);
    
        
        figure
        plot(Data.Nspace(:,1),Data.Nspace(:,2),'o'); hold on
        plot(Data.Nspace(ixSignal_Final,1),Data.Nspace(ixSignal_Final,2),'o','MarkerFaceColor',[0 0 0])
        
        % construct new null model
        P = Data.ExpA(ixSignal_Final,ixSignal_Final); % extract relevant part of null model
        B = Asignal_final - P;          % create a signal modularity matrix
        C = kmeansSweep(Data.Nspace(ixSignal_Final,:),NNegDims(iF)+1,NNegDims(iF)+1,nreps,dims);  % find groups in embedding dimensions: sweep from L to M
        m = sum(sum(Asignal_final))/2;    % number of unique links (or total unique weights)

        for iQ = 1:size(C,2)
            Q(iQ) = computeQ(C(:,iQ),B,m); % compute modularity Q for each clustering
        end
        ixQ = find(Q == min(Q));
        
        Grp_Neg = C(:,ixQ(1));
        [H,C,I] = plotClusterMap(Asignal,Grp_Neg,[],[],'S'); 
        plotorder = Rejection.ixSignal(I);
        title(fnames{iF})
        exportPPTfig(gcf,[fnames{iF} '_QmaxClusterMap'],[10 15 12 12])
    end
end