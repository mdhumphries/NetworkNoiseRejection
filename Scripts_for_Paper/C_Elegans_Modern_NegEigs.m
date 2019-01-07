%% quick script to look at modern C Elegans network's negative eigenvalue results
clearvars; close all

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

%% load rejected data
fname = 'cElegAdjMatAllSynapUndirected';

load(['../Results/Rejected_' fname])

%% visualise the eigenvector
[Nvec,ix] = sort(Data.Nspace(:,1));

% ytick = find(abs(Nvec) > prctile(abs(Nvec),90));
ytick = [1:3 numel(Nvec)-3:numel(Nvec)];
figure
barh(Nvec)
ylabel('Weight')
set(gca,'YTick',ytick,'YTickLabel',Data.nodelabels(ix(ytick),:))
    
%% load null model eigenvectors
% these data were saved from ugly hack of
% "Save_Vectors_All_Data_Networks_Script.m"

load([storepath '/NullModel_Eigenspectrum_' fname])

%% do node rejection in this negative-D space
B = Data.A - Data.ExpA;  % modularity matrix using chosen null model
% NB parameters have been loaded alongside Data

% use rejection parts, but on lower-bound...
optionsReject.Bounds = 'Lower';
Results.Rejection = NodeRejection(B,Data.Emodel,pars.I,Vmodel,optionsReject); % N.B. also calls LowDSpace function to find projections   
Results.NetN = numel(Nvec);

% proportions rejected
Results.PropNodesRetained = numel(Rejection.ixSignal) ./ numel(Nvec);

%% make into clusters
Asignal = Data.A(Results.Rejection.ixSignal,Results.Rejection.ixSignal); 


Results.Asignal = Asignal;
% if 1 eigenvector, just +/- entries after node rejection
Results.Grp_Neg = ones(numel(Results.Rejection.ixSignal),1);
Results.Grp_Neg(Data.Nspace(Results.Rejection.ixSignal) > 0) = 2;

[H,C,I] = plotClusterMap(Results.Asignal,Results.Grp_Neg,[],[],'S'); 
 plotorder = Results.Rejection.ixSignal(I);

% Add node labels
set(gca,'Ytick',1:length(Results.Rejection.ixSignal));
set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
% set(gca,'XTickLabelRotation',90);

