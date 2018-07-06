% script to visualise synthetic model results
% Mark Humphries 8/8/2017

clear all; close all;

addpath('../Helper_Functions/')

fname = 'P_rejection_SyntheticEqual_Noise_20180611T132723';  % P(between) = 0.05

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

nP = numel(Model.P_of_noise);    % P noise
nBatch = size(Results.Time,2);   % number of model repeats
nF = numel(Model.F_noise);

% plotting parameters
I = 0.99; % 99% CI

% colormaps
perf_cmap = brewermap(50,'Blues');
VI_cmap = brewermap(10,'Reds');
colors.line = [0.8 0.7 0.7];
colors.error = [0.8 0.7 0.7];
colors.WCM = [0.8 0.5 0.5];
colors.Config= [0.5 0.5 0.8];
colors.WCMpos = [0.9 0.7 0.7];
colors.Configpos= [0.7 0.7 0.9];
colors.truth = [0.7 0.7 0.7];

M = 5;

%% proportion modular...
figure
[h,xplot,yplot,xtick,ytick] = plotMatrix(Model.F_noise,Model.P_of_noise,Results.ProportionModular.SpectraWCM ./nBatch,[0 0]);
colormap(perf_cmap);
drawgrid(gca,Model.F_noise,Model.P_of_noise,[xplot(1) xplot(end)],[yplot(1) yplot(end)],[0.7 0.7 0.7],1)
xlabel('Size of noise halo')
ylabel('P(noise)')
text(xplot(end),ytick(1),'P<P(out)')
text(xplot(end),ytick(2),'P(out)=P<P(in)')
text(xplot(end),ytick(3),'P(out)<P<P(in)')
text(xplot(end),ytick(4),'P = P(in)')
text(xplot(end),ytick(5),'P(in)<P')
title('Proportion modular')

figure
line([Model.P_of_noise(4) Model.P_of_noise(4)],[0 1],'Color',colors.error); hold on
for iF = 1:numel(Model.F_noise)
    plot(Model.P_of_noise,Results.ProportionModular.SpectraWCM(:,iF) ./nBatch,...
        'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); 
end
ylabel('Proportion modular')
xlabel('P(noise)')
text(Model.P_of_noise(4)-0.025,1.05,'P(noise)=P(in)')

%% number of groups
mWCM = mean(Results.SpectraWCM.Groups,3);

[bndsWCM.L,bndsWCM.U] = bounds(Results.SpectraWCM.Groups,2);
[bndsPosWCM.L,bndsPosWCM.U] = bounds(Results.PosEigWCM.Groups,2);


figure
[h,xplot,yplot,xtick,ytick] = plotMatrix(Model.F_noise,Model.P_of_noise,mWCM,[0 0]);
colormap(perf_cmap);
drawgrid(gca,Model.F_noise,Model.P_of_noise,[xplot(1) xplot(end)],[yplot(1) yplot(end)],[0.7 0.7 0.7],1)
xlabel('Size of noise halo')
ylabel('P(noise)')
text(xplot(end),ytick(1),'P<P(out)')
text(xplot(end),ytick(2),'P(out)=P<P(in)')
text(xplot(end),ytick(3),'P(out)<P<P(in)')
text(xplot(end),ytick(4),'P = P(in)')
text(xplot(end),ytick(5),'P(in)<P')

figure
line([Model.P_of_noise(4) Model.P_of_noise(4)],[0 max(max(mWCM))],'Color',colors.error); hold on
for iF = 1:numel(Model.F_noise)
    plot(Model.P_of_noise,mWCM(:,iF),...
        'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); 
end
ylabel('Number of groups')
xlabel('P(noise)')
text(Model.P_of_noise(4)-0.025,max(max(mWCM))+0.05,'P(noise)=P(in)')

%% accuracy of rejection of nodes
mSpec = mean(Results.SpectraWCM.Specificity,3);
mSens = mean(Results.SpectraWCM.Sensitivity,3);
figure
[h,xplot,yplot,xtick,ytick] = plotMatrix(Model.F_noise,Model.P_of_noise,mSens,[0 0]);
colormap(perf_cmap);
drawgrid(gca,Model.F_noise,Model.P_of_noise,[xplot(1) xplot(end)],[yplot(1) yplot(end)],[0.7 0.7 0.7],1)
xlabel('Size of noise halo')
ylabel('P(noise)')
title('Sensitivity (TPR): TP / TargetP ')
text(xplot(end),ytick(1),'P<P(out)')
text(xplot(end),ytick(2),'P(out)=P<P(in)')
text(xplot(end),ytick(3),'P(out)<P<P(in)')
text(xplot(end),ytick(4),'P = P(in)')
text(xplot(end),ytick(5),'P(in)<P')

figure
line([Model.P_of_noise(4) Model.P_of_noise(4)],[0 1],'Color',colors.error); hold on
for iF = 1:numel(Model.F_noise)
    plot(Model.P_of_noise,mSens(:,iF),...
        'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); 
end
ylabel('Sensitivity: : TP / TargetP')
xlabel('P(noise)')
text(Model.P_of_noise(4)-0.025,1.05,'P(noise)=P(in)')



figure
[h,xplot,yplot,xtick,ytick] = plotMatrix(Model.F_noise,Model.P_of_noise,mSpec,[0 0]);
colormap(perf_cmap);
drawgrid(gca,Model.F_noise,Model.P_of_noise,[xplot(1) xplot(end)],[yplot(1) yplot(end)],[0.7 0.7 0.7],1)

xlabel('Size of noise halo')
ylabel('P(noise)')
title('Specificity (TNR): TN / TargetN')
text(xplot(end),ytick(1),'P<P(out)')
text(xplot(end),ytick(2),'P(out)=P<P(in)')
text(xplot(end),ytick(3),'P(out)<P<P(in)')
text(xplot(end),ytick(4),'P = P(in)')
text(xplot(end),ytick(5),'P(in)<P')

figure
line([Model.P_of_noise(4) Model.P_of_noise(4)],[0 1],'Color',colors.error); hold on
for iF = 1:numel(Model.F_noise)
    plot(Model.P_of_noise,mSpec(:,iF),...
        'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); 
end
ylabel('Specificity (TNR): TN / TargetN')
xlabel('P(noise)')
text(Model.P_of_noise(4)-0.025,1.05,'P(noise)=P(in)')

%% VI of clustering
if blnCluster
    % load clustering
    load([fpath 'Clustering' fname]); 
    
    nBatch = size(ClustResults,1);

    % convert into arrays...   
    for iB = 1:nBatch
        for iP = 1:numel(Model.P_of_noise)
            for iF = 1:numel(Model.F_noise)
                % groups
                nGrps.Louvain(iP,iF,iB) = ClustResults(iB).nGrpsLouvain(iP,iF);
                nGrps.Multiway(iP,iF,iB) = ClustResults(iB).nGrpsMultiway(iP,iF);
                nGrps.QmaxSpectra(iP,iF,iB) = ClustResults(iB).nGrpsQmaxSpectra(iP,iF);
                nGrps.ConsensusSpectra(iP,iF,iB) = ClustResults(iB).nGrpsConsensusSpectra(iP,iF); 
                % VI - noise nodes in solo groups
                VI.Solo.Louvain(iP,iF,iB) = ClustResults(iB).normVILouvainOwn(iP,iF);
                VI.Solo.Multiway(iP,iF,iB) = ClustResults(iB).normVIMultiwayOwn(iP,iF);
                VI.Solo.QmaxSpectra(iP,iF,iB) = ClustResults(iB).normVIQmaxSpectraOwn(iP,iF);
                VI.Solo.ConsensusSpectra(iP,iF,iB) = ClustResults(iB).normVIConsensusSpectraOwn(iP,iF);
                % VI - noise nodes in a single group
                VI.One.Louvain(iP,iF,iB) = ClustResults(iB).normVILouvainOne(iP,iF);
                VI.One.Multiway(iP,iF,iB) = ClustResults(iB).normVIMultiwayOne(iP,iF);
                VI.One.QmaxSpectra(iP,iF,iB) = ClustResults(iB).normVIQmaxSpectraOne(iP,iF);
                VI.One.ConsensusSpectra(iP,iF,iB) = ClustResults(iB).normVIConsensusSpectraOne(iP,iF);            
            end
        end
    end
    % proportion of networks with modules?
    
    
    % numbers of modules
    mLouvain = mean(nGrps.Louvain,3);
   [bnds.Louv.L,bnds.Louv.U] = bounds(nGrps.Louvain,3);
   mMulti = mean(nGrps.Multiway,3);
   [bnds.Multi.L,bnds.Multi.U] = bounds(nGrps.Multiway,3);
   mQmax = mean(nGrps.QmaxSpectra,3);
   [bnds.Qmax.L,bnds.Qmax.U] = bounds(nGrps.QmaxSpectra,3);

    figure
    line([0 max(Model.P_of_noise)],[numel(Model.N), numel(Model.N)],'Color',colors.truth); hold on
    plot(Model.P_of_noise,mLouvain,'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); hold on
    for iF = 1:numel(Model.F_noise)
        line([Model.P_of_noise; Model.P_of_noise],[bnds.Louv.L(:,iF)';bnds.Louv.U(:,iF)'],'Color',colors.WCM)
    end
    xlabel('P(noise)')
    ylabel('Number of modules')
    box off
    if blnExport exportPPTfig(gcf,'nGrpsLouvain',[10 15 8 6]); end
    
    plot(Model.P_of_noise,mMulti,'o-','MarkerSize',M,'Color',colors.Config,'MarkerFaceColor',colors.Config);
    for iF = 1:numel(Model.F_noise)
        line([Model.P_of_noise; Model.P_of_noise],[bnds.Multi.L(:,iF)';bnds.Multi.U(:,iF)'],'Color',colors.Config)
    end
    if blnExport exportPPTfig(gcf,'nGrps_Louv_Multiway',[10 15 8 6]); end

    plot(Model.P_of_noise,mQmax,'o-','MarkerSize',M,'Color',colors.line,'MarkerFaceColor',colors.line);
    for iF = 1:numel(Model.F_noise)
        line([Model.P_of_noise; Model.P_of_noise],[bnds.Qmax.L(:,iF)';bnds.Qmax.U(:,iF)'],'Color',colors.line)
    end
    if blnExport exportPPTfig(gcf,'nGrps_Louv_Multiway_Qmax',[10 15 8 6]); end

    % VI with ground truth: solo groups
    mLouvain = mean(VI.Solo.Louvain,3);
   [bnds.Louv.L,bnds.Louv.U] = bounds(VI.Solo.Louvain,3);
    mMulti = mean(VI.Solo.Multiway,3);
   [bnds.Multi.L,bnds.Multi.U] = bounds(VI.Solo.Multiway,3);
    mQmax = mean(VI.Solo.QmaxSpectra,3);
   [bnds.Qmax.L,bnds.Qmax.U] = bounds(VI.Solo.QmaxSpectra,3);

    figure
    plot(Model.P_of_noise,mLouvain,'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); hold on
    for iF = 1:numel(Model.F_noise)
        line([Model.P_of_noise; Model.P_of_noise],[bnds.Louv.L(:,iF)';bnds.Louv.U(:,iF)'],'Color',colors.WCM)
    end
    xlabel('P(noise)')
    ylabel('VI (normalised): Solo')
    set(gca,'YLim',[0 1]);
    if blnExport exportPPTfig(gcf,'VILouvain',[10 15 8 6]); end
    
    plot(Model.P_of_noise,mMulti,'o-','MarkerSize',M,'Color',colors.Config,'MarkerFaceColor',colors.Config);
    for iF = 1:numel(Model.F_noise)
        line([Model.P_of_noise; Model.P_of_noise],[bnds.Multi.L(:,iF)';bnds.Multi.U(:,iF)'],'Color',colors.Config)
    end
    if blnExport exportPPTfig(gcf,'VILouvain_Multiway',[10 15 8 6]); end

    plot(Model.P_of_noise,mQmax,'o-','MarkerSize',M,'Color',colors.line,'MarkerFaceColor',colors.line);
    for iF = 1:numel(Model.F_noise)
        line([Model.P_of_noise; Model.P_of_noise],[bnds.Qmax.L(:,iF)';bnds.Qmax.U(:,iF)'],'Color',colors.line)
    end
    if blnExport exportPPTfig(gcf,'VILouvain_Multiway_Qmax',[10 15 8 6]); end
    
    % VI with ground truth: one group
    mLouvain = mean(VI.One.Louvain,3);
   [bnds.Louv.L,bnds.Louv.U] = bounds(VI.One.Louvain,3);
    mMulti = mean(VI.One.Multiway,3);
   [bnds.Multi.L,bnds.Multi.U] = bounds(VI.One.Multiway,3);
    mQmax = mean(VI.One.QmaxSpectra,3);
   [bnds.Qmax.L,bnds.Qmax.U] = bounds(VI.One.QmaxSpectra,3);

    figure
    plot(Model.P_of_noise,mLouvain,'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); hold on
    for iF = 1:numel(Model.F_noise)
        line([Model.P_of_noise; Model.P_of_noise],[bnds.Louv.L(:,iF)';bnds.Louv.U(:,iF)'],'Color',colors.WCM)
    end
    xlabel('P(noise)')
    ylabel('VI (normalised): One')
    set(gca,'YLim',[0 1]);
    if blnExport exportPPTfig(gcf,'VILouvain',[10 15 8 6]); end
    
    plot(Model.P_of_noise,mMulti,'o-','MarkerSize',M,'Color',colors.Config,'MarkerFaceColor',colors.Config);
    for iF = 1:numel(Model.F_noise)
        line([Model.P_of_noise; Model.P_of_noise],[bnds.Multi.L(:,iF)';bnds.Multi.U(:,iF)'],'Color',colors.Config)
    end
    if blnExport exportPPTfig(gcf,'VILouvain_Multiway',[10 15 8 6]); end

    plot(Model.P_of_noise,mQmax,'o-','MarkerSize',M,'Color',colors.line,'MarkerFaceColor',colors.line);
    for iF = 1:numel(Model.F_noise)
        line([Model.P_of_noise; Model.P_of_noise],[bnds.Qmax.L(:,iF)';bnds.Qmax.U(:,iF)'],'Color',colors.line)
    end
    if blnExport exportPPTfig(gcf,'VILouvain_Multiway_Qmax',[10 15 8 6]); end
    
end