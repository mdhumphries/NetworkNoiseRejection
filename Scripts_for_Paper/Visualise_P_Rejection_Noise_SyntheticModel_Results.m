% script to visualise synthetic model results
% Mark Humphries 8/8/2017

clear all; close all;

addpath('../Helper_Functions/')

fname = 'P_rejection_SyntheticEqual_Noise_20180607T132417';
fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';

blnCluster = 0;  % if done clustering, set this to 1
blnExport = 0;

% load results and set maps
load([fpath fname]);  

nP = numel(Model.P_of_noise);    % P noise
nBatch = size(Results.Time,2);   % number of model repeats
nF = numel(Model.F_noise);

% plotting parameters
I = 0.99; % 99% CI

% colormaps
perf_cmap = brewermap(10,'Blues');
VI_cmap = brewermap(10,'Reds');
colors.line = [0.8 0.7 0.7];
colors.error = [0.8 0.7 0.7];
colors.WCM = [0.8 0.5 0.5];
colors.Config= [0.5 0.5 0.8];
colors.WCMpos = [0.9 0.7 0.7];
colors.Configpos= [0.7 0.7 0.9];
colors.truth = [0.7 0.7 0.7];

M = 5;

%% number of groups
mWCM = mean(Results.SpectraWCM.Groups,3);

[bndsWCM.L,bndsWCM.U] = bounds(Results.SpectraWCMGroups,2);
[bndsPosWCM.L,bndsPosWCM.U] = bounds(Results.PosEigWCMGroups,2);

figure
imagesc(Model.F_noise,Model.P_of_noise,mWCM)
colormap(perf_cmap);
xlabel('Size of noise halo')
ylabel('P(noise)')
set(gca,'YDir','normal','YTick',Model.P_of_noise,'XTick',Model.F_noise)
text(1.2,Model.P_of_noise(1),'P<P(out)')
text(1.2,Model.P_of_noise(2),'P(out)<P<P(in)')
text(1.2,Model.P_of_noise(3),'P(in)<P')

%% accuracy of rejection of nodes
mSpec = mean(Results.SpectraWCM.Specificity,3);
mSens = mean(Results.SpectraWCM.Sensitivity,3);
figure
imagesc(Model.F_noise,Model.P_of_noise,mSpec)
colormap(perf_cmap);
xlabel('Size of noise halo')
ylabel('P(noise)')
set(gca,'YDir','normal','YTick',Model.P_of_noise,'XTick',Model.F_noise)
title('Sensitivity (TPR): TP / TargetP ')
text(1.2,Model.P_of_noise(1),'P<P(out)')
text(1.2,Model.P_of_noise(2),'P(out)<P<P(in)')
text(1.2,Model.P_of_noise(3),'P(in)<P')

figure
imagesc(Model.F_noise,Model.P_of_noise,mSens)
colormap(perf_cmap);
xlabel('Size of noise halo')
ylabel('P(noise)')
set(gca,'YDir','normal','YTick',Model.P_of_noise,'XTick',Model.F_noise)
title('Specificity (TNR): TN / TargetN')
text(1.2,Model.P_of_noise(1),'P<P(out)')
text(1.2,Model.P_of_noise(2),'P(out)<P<P(in)')
text(1.2,Model.P_of_noise(3),'P(in)<P')

%% VI of clustering
if blnCluster
    % load clustering
    load([fpath 'Clustering' fname]); 
    
    nBatch = size(ClustResults,1);

    % convert into arrays...   
    for iB = 1:nBatch
        for iP = 1:numel(Model.P_of_within)
            % groups
            nGrps.Louvain(iP,iB) = ClustResults(iB).nGrpsLouvain(iP);
            nGrps.Multiway(iP,iB) = ClustResults(iB).nGrpsMultiway(iP);
            nGrps.QmaxSpectra(iP,iB) = ClustResults(iB).nGrpsQmaxSpectra(iP);
            nGrps.ConsensusSpectra(iP,iB) = ClustResults(iB).nGrpsConsensusSpectra(iP); 
            % VI
            VI.Louvain(iP,iB) = ClustResults(iB).normVILouvain(iP);
            VI.Multiway(iP,iB) = ClustResults(iB).normVIMultiway(iP);
            VI.QmaxSpectra(iP,iB) = ClustResults(iB).normVIQmaxSpectra(iP);
            VI.ConsensusSpectra(iP,iB) = ClustResults(iB).normVIConsensusSpectra(iP);             
        end
    end
    % proportion of networks with modules?
    
    
    % numbers of modules
    mLouvain = mean(nGrps.Louvain,2);
   [bnds.Louv.L,bnds.Louv.U] = bounds(nGrps.Louvain,2);
   mMulti = mean(nGrps.Multiway,2);
   [bnds.Multi.L,bnds.Multi.U] = bounds(nGrps.Multiway,2);
    mQmax = mean(nGrps.QmaxSpectra,2);
   [bnds.Qmax.L,bnds.Qmax.U] = bounds(nGrps.QmaxSpectra,2);

    figure
    line([0 max(Pdiff)],[numel(Model.N), numel(Model.N)],'Color',colors.truth); hold on
    plot(Pdiff,mLouvain,'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); hold on
    line([Pdiff; Pdiff],[bnds.Louv.L';bnds.Louv.U'],'Color',colors.WCM)
    xlabel('P(within) - P(between)')
    ylabel('Number of modules')
    box off
    if blnExport exportPPTfig(gcf,'nGrpsLouvain',[10 15 8 6]); end
    
    plot(Pdiff,mMulti,'o-','MarkerSize',M,'Color',colors.Config,'MarkerFaceColor',colors.Config);
    line([Pdiff; Pdiff],[bnds.Multi.L';bnds.Multi.U'],'Color',colors.Config)
    if blnExport exportPPTfig(gcf,'nGrps_Louv_Multiway',[10 15 8 6]); end

    plot(Pdiff,mQmax,'o-','MarkerSize',M,'Color',colors.line,'MarkerFaceColor',colors.line);
    line([Pdiff; Pdiff],[bnds.Qmax.L';bnds.Qmax.U'],'Color',colors.line)
    if blnExport exportPPTfig(gcf,'nGrps_Louv_Multiway_Qmax',[10 15 8 6]); end

    % VI with ground truth
    mLouvain = mean(VI.Louvain,2);
   [bnds.Louv.L,bnds.Louv.U] = bounds(VI.Louvain,2);
   mMulti = mean(VI.Multiway,2);
   [bnds.Multi.L,bnds.Multi.U] = bounds(VI.Multiway,2);
    mQmax = mean(VI.QmaxSpectra,2);
   [bnds.Qmax.L,bnds.Qmax.U] = bounds(VI.QmaxSpectra,2);

    figure
    plot(Pdiff,mLouvain,'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); hold on
    line([Pdiff; Pdiff],[bnds.Louv.L';bnds.Louv.U'],'Color',colors.WCM)
    xlabel('P(within) - P(between)')
    ylabel('VI (normalised)')
    set(gca,'YLim',[0 1]);
    if blnExport exportPPTfig(gcf,'VILouvain',[10 15 8 6]); end
    
    plot(Pdiff,mMulti,'o-','MarkerSize',M,'Color',colors.Config,'MarkerFaceColor',colors.Config);
    line([Pdiff; Pdiff],[bnds.Multi.L';bnds.Multi.U'],'Color',colors.Config)
    if blnExport exportPPTfig(gcf,'VILouvain_Multiway',[10 15 8 6]); end

    plot(Pdiff,mQmax,'o-','MarkerSize',M,'Color',colors.line,'MarkerFaceColor',colors.line);
    line([Pdiff; Pdiff],[bnds.Qmax.L';bnds.Qmax.U'],'Color',colors.line)
    if blnExport exportPPTfig(gcf,'VILouvain_Multiway_Qmax',[10 15 8 6]); end

end