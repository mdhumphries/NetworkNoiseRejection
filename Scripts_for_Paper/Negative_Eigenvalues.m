%% negative eigenvalues
%
% Mark Humphries 26/7/2018

clearvars; close all

load('../Results/Network_Rejection_Table')

% location of null model eigenvectors
storepath = 'C:\Users\lpzmdh\Dropbox\Analyses\Networks\DataNets_Null_EigVectors\';

addpath('../Network_Spectra_Functions/')
addpath('../Network_Analysis_Functions/')

%% find networks with negative eigenvalues below lower limit
fnames = Network_Rejection_Table.NetworkName(Network_Rejection_Table.SparseWCM_NegDn > 0);
NNegDims = Network_Rejection_Table.SparseWCM_NegDn(Network_Rejection_Table.SparseWCM_NegDn > 0);

for iF = 1:numel(fnames)
    % get rejection results
    load(['../Results/Rejected_' fnames{iF}])
    
    % plot eigenvector
    [Nvec,ix] = sort(Data.Nspace(:,1));
    
    ytick = find(abs(Nvec) > prctile(abs(Nvec),90));
    
    figure
    barh(Nvec)
    ylabel('Weight')
    set(gca,'YTick',ytick,'YTickLabel',Data.nodelabels(ix(ytick),:))
    title(fnames{iF})
    
    if NNegDims(iF) > 1
        [Nvec,ix] = sort(Data.Nspace(:,2));
        ytick = find(abs(Nvec) > prctile(abs(Nvec),90));
        figure
        barh(Nvec)
        ylabel('Weight')
        set(gca,'YTick',ytick,'YTickLabel',Data.nodelabels(ix(ytick),:))
        title([fnames{iF} ' 2nd negative eigenvector'])   
    end
    
    %% load null model eigenvectors
    load([storepath '/NullModel_Eigenspectrum_' fnames{iF}])
    
    %% do node rejection in this negative-D space
    B = Data.A - Data.ExpA;  % modularity matrix using chosen null model
    % NB parameters have been loaded alongside Data
    
    % use rejection parts, but on lower-bound...
    optionsReject.Bounds = 'Lower';
    Rejection = NodeRejection(B,Data.Emodel,pars.I,Vmodel,optionsReject); % N.B. also calls LowDSpace function to find projections

    keyboard
    
    
    
    %% make into clusters
    
    % if 1 eigenvector, just +/- entries after node rejection
    Grp_Neg = ones(numel(Rejection.ixSignal),1);
    Grp_Neg(Data.Nspace(Rejection.ixSignal) > 0) = 2;
    Asignal = Data.A(Rejection.ixSignal,Rejection.ixSignal); 
    
    %% if 2 or more, pass to Qmax....
    
    % refine signal matrix for clustering
%     % connected signal matrix: find largest component, and use that - store
%     % others
%     [Asignal_comp,ixRetain,~,~] = prep_A(Asignal); 
%     ixSignal_comp = Rejection.ixSignal(ixRetain);  % original node indices
% 
%     % and then strip out leaves - nodes with single links
%     K = sum(Asignal_comp);
%     ixLeaves = find(K==1); ixKeep = find(K > 1);
% 
%     ixSignal_Final = ixSignal_comp(ixKeep);
%     ixSignal_Leaves = ixSignal_comp(ixLeaves);
%     Asignal_final = Asignal_comp(ixKeep,ixKeep);

end