%% script to visualise Full and Signal networks
clear all; close all

fname = 'Lesmis.mat'; 
blnVizNet = 1;  % network visualisation - if MATLAB BGL installed, appropriate for platform:
fontsize = 6;

if blnVizNet
    % Traud Mucha Porter visualisation tools
    addpath('Traud_Mucha_Porter_CommunityVisualisation/');

    % needs MATLAB BGL Toolbox on your path - change to your local path
    % here:
    if ismac
        bglpath = genpath('/Users/mqbssmhg/Dropbox/My Toolboxes/Graph_theory/matlab_bglOSX64/');  % generate path to local BGL and all its subdirectories
    else
        bglpath = genpath('C:\Users\mqbssmhg.DS\Dropbox\My Toolboxes\Graph_theory\matlab_bgl\');
    end
    % add to current MATLAB path
    addpath(bglpath); 
end

%% Image plot of A, with labels
k = sum(A);
[~,iK] = sort(k,'descend'); 
cmap = brewermap(10,'Greys');
figure; 
imagesc(A(iK,iK)); colormap(cmap);
set(gca,'Ytick',1:length(A));
set(gca,'Yticklabel',nodelabels(iK),'Fontsize',fontsize);
