%% negative eigenvalues
%
% Mark Humphries 26/7/2018

clearvars; close all

load('../Results/Network_Rejection_Table')

% find networks with negative eigenvalues below lower limit
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
    
    %% make into two clusters
    
    
    
end