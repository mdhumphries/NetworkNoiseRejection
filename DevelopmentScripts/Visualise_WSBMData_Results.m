%% script to visualise the synthetic test data
% 3 variables:
% 1. proportion of network that is not in communities
% 2. proportion of within weight used between communities
% 3. proportion of between weights used from periphery to communities
%
% Mark Humphries 21/7/2017

close all; clear all;

fname = 'WSBM_20170721T170120';

load(['../Results/' fname]);  % Main table is 'core_perf'

nNoise= numel(Communities.cp_gg_grid);  % from periphery to core communities
nMix= numel(Communities.gg_cc_grid);    % between communities

% for 3 fractions of network (lowest, middle, highest)

% colormaps
perf_cmap = brewermap(10,'*Blues');
VI_cmap = brewermap(10,'Reds');

%% for each fraction of network
for iP = Communities.fraction_periphery_grid
    ixTable = core_perf.fraction_periphery == iP;
    
    % plot performance as function of mixing and noise [2D heatmap]: 1 is best
    matSignal= reshape(core_perf.signal(ixTable),nNoise,nMix);
    matAdditional= reshape(core_perf.additional(ixTable),nNoise,nMix);
    matIndexMix = reshape(core_perf.mixing_index(ixTable),nNoise,nMix);
    matIndexNoise = reshape(core_perf.noise_index(ixTable),nNoise,nMix);

    figure
    subplot(121)
    colormap(perf_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matSignal);
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Recovery of signal. Periphery = ' num2str(iP)])

    subplot(122)
    colormap(perf_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matAdditional);
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Additional periphery. Periphery = ' num2str(iP)])
    colorbar
    
    %% plot VI as a function of mixing and noise [2D heatmap]: 0 is best
    matVIQmax = reshape(core_perf.nVIQmaxSpectra(ixTable),nNoise,nMix);
    matVICons = reshape(core_perf.nVIConsensusSpectra(ixTable),nNoise,nMix);
    matVILouvainBest = reshape(core_perf.nVILouvainMin(ixTable),nNoise,nMix);
    matVIMultiway = reshape(core_perf.nVIMultiway(ixTable),nNoise,nMix);

    matFullVIQmax = reshape(core_perf.nVIFull_QmaxSpectra(ixTable),nNoise,nMix);
    matFullVICons = reshape(core_perf.nVIFull_ConsensusSpectra(ixTable),nNoise,nMix);
    matFullVILouvainBest = reshape(core_perf.nVIFull_LouvainMin(ixTable),nNoise,nMix);
    matFullVIMultiway = reshape(core_perf.nVIFull_Multiway(ixTable),nNoise,nMix);

    
    % one per clustering method (Qmax, Consensus, best Louvain, multi-way
    % vector)
    
    figure
    % Q max
    subplot(421)
    colormap(VI_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matVIQmax);
    colorbar
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Q_max VI on Signal network. Periphery = ' num2str(iP)])
    
    subplot(422)
    colormap(VI_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matFullVIQmax);
    colorbar
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Q_max VI on Full network. Periphery = ' num2str(iP)])
    
    
    % Consensus
    subplot(423)
    colormap(VI_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matVICons);
    colorbar
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Consensus VI on Signal network. Periphery = ' num2str(iP)])
    
    subplot(424)
    colormap(VI_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matFullVICons);
    colorbar
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Consensus VI on Full network. Periphery = ' num2str(iP)])

    
    % Best Louvain
    subplot(425)
    colormap(VI_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid, matVILouvainBest);
    colorbar
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Best Louvain VI on Signal network. Periphery = ' num2str(iP)])
    
    subplot(426)
    colormap(VI_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matFullVILouvainBest);
    colorbar
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Best Louvain VI on Full network. Periphery = ' num2str(iP)])

     % Multi-way
    subplot(427)
    colormap(VI_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matVIMultiway);
    colorbar
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Multiway VI on Signal network. Periphery = ' num2str(iP)])
    
    subplot(428)
    colormap(VI_cmap)
    imagesc(Communities.gg_cc_grid,Communities.cp_gg_grid,matFullVIMultiway );
    colorbar
    set(gca,'YDir','normal')
    ylabel('Noise')
    xlabel('Mixing')
    title(['Multiway VI on Full network. Periphery = ' num2str(iP)])
   
    h=findobj(gcf,'type','axes');
    q=get(h,'clim');
    set(h,'clim',[min([q{:}]) max([q{:}])])

end