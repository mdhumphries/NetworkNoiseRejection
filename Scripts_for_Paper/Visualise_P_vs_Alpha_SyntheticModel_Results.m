%% script to visualise synthetic model results
% Mark Humphries 8/8/2017

clear all; close all;

addpath('../Helper_Functions/')

fname = 'SyntheticEqual_NoNoise_20180604T150903';
blnCluster = 0;  % if done clustering, set this to 1

% load results and set maps
load(['../Results/' fname]);  % Main table is 'core_perf'

nAlpha = numel(Model.alpha_range);  % split of weights between within and between
nP = numel(Model.P_of_within);    % within communities P
nBatch = size(Results.Time,3);   % number of model repeats

% colormaps
perf_cmap = brewermap(10,'Blues');
VI_cmap = brewermap(10,'Reds');

%% rescaling of axes for Pcolor plots
xsc = [Model.alpha_range 2*Model.alpha_range(end) - Model.alpha_range(end-1)];
xtick = Model.alpha_range + diff(xsc)/2;

ysc = [Model.P_of_within 2*Model.P_of_within(end) - Model.P_of_within(end-1)];
ytick = Model.P_of_within + diff(ysc)/2;

% this part of rescaling solution from: https://uk.mathworks.com/matlabcentral/newsreader/view_thread/30404
% xsc = (xsc-Model.alpha_range(1))/(xsc(1)-xsc(end))*(Model.alpha_range(1)-Model.alpha_range(end)) + Model.alpha_range(1); % scale proportionally
% ysc = (ysc-Model.P_of_within(1))/(ysc(1)-ysc(end))*(Model.P_of_within(1)-Model.P_of_within(end)) + Model.P_of_within(1); % scale proportionally


%% processing time
meanTime = mean(Results.Time,3);
% padding for PColor function
plotmap = [meanTime meanTime(:,end)]; 
plotmap = [plotmap; plotmap(end,:)];

figure
colormap(perf_cmap)
pcolor(xsc,ysc,plotmap);
set(gca,'XTick',xtick,'XTickLabel',Model.alpha_range)
set(gca,'YTick',ytick,'YTickLabel',Model.P_of_within)
% pcolor(Model.alpha_range,Model.P_of_within,Results.Time);
colorbar
xlabel('Proportion of strength within community')
ylabel('P(within)')
title('Time to do all processing (s)')

%% recovered number of groups 
figure

subplot(221),
plotmap = [Results.SpectraWCMGroups Results.SpectraWCMGroups(:,end)]; 
plotmap = [plotmap; plotmap(end,:)];
colormap(perf_cmap)
pcolor(xsc,ysc,plotmap);
set(gca,'XTick',xtick,'XTickLabel',Model.alpha_range)
set(gca,'YTick',ytick,'YTickLabel',Model.P_of_within)
colorbar
xlabel('Proportion of strength within community')
ylabel('P(within)')
title('Number of groups (WCM Spectra)')

subplot(222),
plotmap = [Results.SpectraConfigGroups Results.SpectraConfigGroups(:,end)]; 
plotmap = [plotmap; plotmap(end,:)];
colormap(perf_cmap)
pcolor(xsc,ysc,plotmap);
set(gca,'XTick',xtick,'XTickLabel',Model.alpha_range)
set(gca,'YTick',ytick,'YTickLabel',Model.P_of_within)
colorbar
xlabel('Proportion of strength within community')
ylabel('P(within)')
title('Number of groups (Config Model Spectra)')

subplot(223),
plotmap = [Results.PosEigWCMGroups Results.PosEigWCMGroups(:,end)]; 
plotmap = [plotmap; plotmap(end,:)];
colormap(perf_cmap)
pcolor(xsc,ysc,plotmap);
set(gca,'XTick',xtick,'XTickLabel',Model.alpha_range)
set(gca,'YTick',ytick,'YTickLabel',Model.P_of_within)
colorbar
xlabel('Proportion of strength within community')
ylabel('P(within)')
title('Number of groups (WCM all positive eigs)')

subplot(224),
plotmap = [Results.PosEigConfigGroups Results.PosEigConfigGroups(:,end)]; 
plotmap = [plotmap; plotmap(end,:)];
colormap(perf_cmap)
pcolor(xsc,ysc,plotmap);
set(gca,'XTick',xtick,'XTickLabel',Model.alpha_range)
set(gca,'YTick',ytick,'YTickLabel',Model.P_of_within)
colorbar
xlabel('Proportion of strength within community')
ylabel('P(within)')
title('Number of groups (Config all positive eigs)')


%% VI of clustering
if blnCluster
    % load clustering
    load(['../Results/Clustering' fname]);  % Main table is 'core_perf'

    % plot clustering
    figure

    subplot(221),
    plotmap = [Results.nVIFull_QmaxSpectra Results.nVIFull_QmaxSpectra(:,end)]; 
    plotmap = [plotmap; plotmap(end,:)];
    colormap(VI_cmap)
    pcolor(xsc,ysc,plotmap);
    set(gca,'XTick',xtick,'XTickLabel',Model.alpha_range)
    set(gca,'YTick',ytick,'YTickLabel',Model.P_of_within)
    colorbar
    xlabel('Proportion of strength within community')
    ylabel('P(within)')
    title('VI Q max solution')

    subplot(222),
    plotmap = [Results.nVIFull_ConsensusSpectra Results.nVIFull_ConsensusSpectra(:,end)]; 
    plotmap = [plotmap; plotmap(end,:)];
    colormap(VI_cmap)
    pcolor(xsc,ysc,plotmap);
    set(gca,'XTick',xtick,'XTickLabel',Model.alpha_range)
    set(gca,'YTick',ytick,'YTickLabel',Model.P_of_within)
    colorbar
    xlabel('Proportion of strength within community')
    ylabel('P(within)')
    title('VI Consensus')

    subplot(223),
    plotmap = [Results.nVIFull_LouvainMin Results.nVIFull_LouvainMin(:,end)]; 
    plotmap = [plotmap; plotmap(end,:)];
    colormap(VI_cmap)
    pcolor(xsc,ysc,plotmap);
    set(gca,'XTick',xtick,'XTickLabel',Model.alpha_range)
    set(gca,'YTick',ytick,'YTickLabel',Model.P_of_within)
    colorbar
    xlabel('Proportion of strength within community')
    ylabel('P(within)')
    title('VI best Louvain')

    subplot(224),
    plotmap = [Results.nVIFull_Multiway Results.nVIFull_Multiway(:,end)]; 
    plotmap = [plotmap; plotmap(end,:)];
    colormap(VI_cmap)
    pcolor(xsc,ysc,plotmap);
    set(gca,'XTick',xtick,'XTickLabel',Model.alpha_range)
    set(gca,'YTick',ytick,'YTickLabel',Model.P_of_within)
    colorbar
    xlabel('Proportion of strength within community')
    ylabel('P(within)')
    title('VI Multiway')

    % recale to all have the same colour-scale
    h=findobj(gcf,'type','axes');
    q=get(h,'clim');
    set(h,'clim',[min([q{:}]) max([q{:}])])
end