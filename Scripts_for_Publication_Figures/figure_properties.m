% figure properties for PLoS Biology paper

format = 'png'; % for panels tricky for EPS (e.g. Pcolor plots)
color = 'rgb';
dpi = 600;
fontsize = 8;
fontname = 'Arial';
M = 4; % marker size for univariate scatter plots
sym = 'o';  % markers for scatters and strip plots

Units = 'centimeters';

% line widths
widths.plot = 0.75;
widths.error = 0.5;
widths.axis = 0.5;


% panel sizes
figsize = [4 4];

% colours for matrix checkerboards (if using)
perf_cmap = brewermap(50,'Blues');
VI_cmap = brewermap(10,'Reds');

% colours for lins
colors.truth = [0 0 0];
colors.selecteigs = [0.8 0.5 0.7];
colors.poseigs = [0.7 0.7 0.7];

colors.Louvain = [0.8 0.5 0.5];
colors.Multiway = [0.5 0.5 0.8];
colors.Qmax = [0.8 0.7 0.7];

% plotting ticks and limits for binsize plots

% exportpath
if ispc
    exportpath = 'C:/Users/lpzmdh/Dropbox/My Papers/Networks/Noise rejection for networks/Figures/Panels/';
else
    exportpath = '/Users/mqbssmhg/Dropbox/My Papers/Networks/Noise rejection for networks/Figures/Panels/';
end