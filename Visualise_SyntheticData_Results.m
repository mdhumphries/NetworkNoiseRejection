%% script to visualise the synthetic test data
% 3 variables:
% 1. proportion of network that is not in communities
% 2. proportion of within weight used between communities
% 3. proportion of between weights used from periphery to communities
%
% Mark Humphries 21/7/2017

close all; clear all;

fname = 'WSBM_20170721T115539';

load(['Results/' fname]);  % Main table is 'core_perf'

nNoise= numel(Communities.cp_gg_grid);  % from periphery to core communities
nMix= numel(Communities.gg_cc_grid);    % between communities

% for 3 fractions of network (lowest, middle, highest)

% colormaps
perf_cmap = brewermap(10,'*Blues');
VI_cmap = brewermap(10,'Reds');

for iP = Communities.fraction_periphery_grid
    ixTable = core_perf.fraction_periphery == iP;
    
    % plot performance as function of mixing and noise [2D heatmap]: 1 is best
    matPerform = reshape(core_perf.performance(ixTable),nNoise,nMix);
    matIndexMix = reshape(core_perf.mixing_index(ixTable),nNoise,nMix);
    matIndexNoise = reshape(core_perf.noise_index(ixTable),nNoise,nMix);

    figure
    colormap(perf_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matPerform);
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Recovery of signal. Periphery = ' num2str(iP)])

    % plot VI as a function of mixing and noise [2D heatmap]: 0 is best
    matVIQmax = reshape(core_perf.nVIQmaxSpectra(ixTable),nNoise,nMix);
    matVICons = reshape(core_perf.nVIConsensusSpectra(ixTable),nNoise,nMix);
    matVILouvainBest = reshape(core_perf.nVILouvainMin(ixTable),nNoise,nMix);
    matVIMultiway = reshape(core_perf.nVIMultiway(ixTable),nNoise,nMix);

    % one per clustering method (Qmax, Consensus, best Louvain, multi-way
    % vector)
    
    figure
    colormap(VI_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matVIQmax);
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Q_max VI. Periphery = ' num2str(iP)])
    
end