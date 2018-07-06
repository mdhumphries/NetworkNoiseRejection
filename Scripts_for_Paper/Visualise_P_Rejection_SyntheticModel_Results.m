%% script to visualise synthetic model results
% Mark Humphries 8/8/2017

clear all; close all;

addpath('../Helper_Functions/')

% fname = 'P_rejection_SyntheticEqual_NoNoise_20180702T124726'; % P(between) = 0.05
% fname = 'P_rejection_SyntheticEqual_NoNoise_20180705T141726'; % P(between) = 0.15
fname = 'P_rejection_SyntheticUnequal_NoNoise_20180705T180300'; % P(between) = 0.05 and 3 unequal groups

% path to storage of huge analysis files
if ispc
    fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';
else
    fpath = '/Users/mqbssmhg/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';
end

blnCluster = 1;  % if done clustering, set this to 1
blnExport = 0;

% load results and set maps
load([fpath fname]);  

nP = numel(Model.P_of_within);    % within communities P
nBatch = size(Results.Time,2);   % number of model repeats

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

Pdiff = Model.P_of_within - Model.P.between;

%% processing time
meanTime = mean(Results.Time,2);
CITime = CIfromSEM(std(Results.Time,[],2),ones(nP,1).*nBatch,I);

figure
plot(Pdiff,meanTime,'o-','Color',colors.line,'MarkerSize',M,'MarkerFaceColor',colors.line)
line([Pdiff; Pdiff],[(meanTime-CITime)';(meanTime+CITime)'],'Color',colors.error)
set(gca,'XTick',Pdiff,'XTickLabel',Pdiff)
ylabel('Time (s)')
xlabel('P(within) - P(between)')
box off
%% recovered proportions of networks with groups

figure
plot(Pdiff,Results.ProportionModular.SpectraWCM ./nBatch,'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); hold on
plot(Pdiff,Results.ProportionModular.SpectraConfig ./nBatch,'o-','MarkerSize',M,'Color',colors.Config,'MarkerFaceColor',colors.Config)
plot(Pdiff,Results.ProportionModular.PosEigWCM ./nBatch,'o-','MarkerSize',M,'Color',colors.WCMpos,'MarkerFaceColor',colors.WCMpos);
set(gca,'XTick',Pdiff,'XTickLabel',Pdiff)
xlabel('P(within) - P(between)')
ylabel('Networks with modules (%)')
box off

%% number of groups
mWCM = mean(Results.SpectraWCMGroups,2);
mConfig = mean(Results.SpectraConfigGroups,2);
mPosWCM = mean(Results.PosEigWCMGroups,2);
[bndsWCM.L,bndsWCM.U] = bounds(Results.SpectraWCMGroups,2);
[bndsConfig.L,bndsConfig.U] = bounds(Results.SpectraConfigGroups,2);
[bndsPosWCM.L,bndsPosWCM.U] = bounds(Results.PosEigWCMGroups,2);

figure
line([0 max(Pdiff)],[numel(Model.N), numel(Model.N)],'Color',colors.truth); hold on
plot(Pdiff,mWCM,'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); hold on
line([Pdiff; Pdiff],[bndsWCM.L';bndsWCM.U'],'Color',colors.WCM)
plot(Pdiff,mConfig,'o-','MarkerSize',M,'Color',colors.Config,'MarkerFaceColor',colors.Config); hold on
line([Pdiff; Pdiff],[bndsConfig.L';bndsConfig.U'],'Color',colors.Config)

plot(Pdiff,mPosWCM,'o-','MarkerSize',M,'Color',colors.WCMpos,'MarkerFaceColor',colors.WCMpos);
line([Pdiff; Pdiff],[bndsPosWCM.L';bndsPosWCM.U'],'Color',colors.WCMpos)

set(gca,'XTick',Pdiff,'XTickLabel',Pdiff)
xlabel('P(within) - P(between)')
ylabel('Number of modules')
box off

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