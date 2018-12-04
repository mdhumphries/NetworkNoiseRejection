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

% hack worm in...
% fnames{end+1} = 'Worm279_Wmatrix';
% NNegDims(end+1) = 1;

for iF = 1:numel(fnames)
    % get rejection results
    load(['../Results/Rejected_' fnames{iF}])
    
    Results(iF).NetName = fnames{iF};
    Results(iF).nodelabels = Data.nodelabels;
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
    Results(iF).Rejection = NodeRejection(B,Data.Emodel,pars.I,Vmodel,optionsReject); % N.B. also calls LowDSpace function to find projections   
    Results(iF).NetN = numel(Nvec);
    
    % proportions rejected
    Results(iF).PropNodesRetained(iF) = numel(Rejection.ixSignal) ./ numel(Nvec);
    
    %% make into clusters
    Asignal = Data.A(Results(iF).Rejection.ixSignal,Results(iF).Rejection.ixSignal); 

    if NNegDims(iF) == 1
        Results(iF).Asignal = Asignal;
        % if 1 eigenvector, just +/- entries after node rejection
        Results(iF).Grp_Neg = ones(numel(Results(iF).Rejection.ixSignal),1);
        Results(iF).Grp_Neg(Data.Nspace(Results(iF).Rejection.ixSignal) > 0) = 2;

        [H,C,I] = plotClusterMap(Results(iF).Asignal,Results(iF).Grp_Neg,[],[],'S'); 
         plotorder = Results(iF).Rejection.ixSignal(I);
    
        % Add node labels
        set(gca,'Ytick',1:length(Results(iF).Rejection.ixSignal));
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
        ixLeaves = find(K==1); 
        ixKeep = find(K > 1);

        Results(iF).ixSignal_Final = ixSignal_comp(ixKeep);
        Results(iF).ixSignal_Leaves = ixSignal_comp(ixLeaves);
        Results(iF).Asignal_final = Asignal_comp(ixKeep,ixKeep);
    
        
        figure
        plot(Data.Nspace(:,1),Data.Nspace(:,2),'o'); hold on
        plot(Data.Nspace(Results(iF).ixSignal_Final,1),Data.Nspace(Results(iF).ixSignal_Final,2),'o','MarkerFaceColor',[0 0 0])
        
        % construct new null model
        P = Data.ExpA(Results(iF).ixSignal_Final,Results(iF).ixSignal_Final); % extract relevant part of null model
        B = Results(iF).Asignal_final - P;          % create a signal modularity matrix
        Results(iF).sweep = kmeansSweep(Data.Nspace(Results(iF).ixSignal_Final,:),NNegDims(iF)+1,NNegDims(iF)+1,nreps,dims);  % find groups in embedding dimensions: sweep from L to M
        m = sum(sum(Results(iF).Asignal_final))/2;    % number of unique links (or total unique weights)

        for iQ = 1:size(Results(iF).sweep,2)
            Results(iF).Q(iQ) = computeQ(Results(iF).sweep(:,iQ),B,m); % compute modularity Q for each clustering
        end
        ixQ = find(Results(iF).Q == min(Results(iF).Q));
        
        Results(iF).Grp_Neg = Results(iF).sweep(:,ixQ(1));
        [H,C,I] = plotClusterMap(Results(iF).Asignal_final,Results(iF).Grp_Neg,[],[],'S'); 
        plotorder = Results(iF).ixSignal_Final(I);
        title(fnames{iF})
        exportPPTfig(gcf,[fnames{iF} '_QmaxClusterMap'],[10 15 12 12])
    end
end

save('../Results/NegEigResults','Results')

%% Pol Blogs: weird stuff!
CCon = makeConsensusMatrix(Results(end).sweep);

uniqueCs = CCon(triu(ones(size(CCon)),1) == 1);
figure   % there are a few entries clustered exactly 50% of the time together
hist(uniqueCs,30);

ids = find(triu(CCon) > 0.4 & triu(CCon) < 0.6);
[rs,cs] = ind2sub(size(CCon),ids);

ixJump = unique(rs);
for i=1:numel(ixJump)
    count(i) = sum(rs == ixJump(i));
end

figure 
bar(count)

ixTri = ixJump(count > median(count));

% make groups
[blnC,G] = CheckConvergenceConsensus(CCon);

% Results(iF).Grp_ConNeg = G;
% [H,C,I] = plotClusterMap(Results(iF).Asignal_final,Results(iF).Grp_ConNeg,[],[],'S'); 
% 
% save('../Results/NegEigResults','Results')









