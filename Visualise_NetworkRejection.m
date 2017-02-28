%% script to visualise Full and Signal networks
% Visualisations of network layouts need:
% (1) Traud-Mucha-Porter toolbox (included in GitHub)
% (2) MATLAB BGL Toolbox:
%           Win32, Win64, Mac32, Linux: https://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/ 
%           Mac64: http://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/old/matlab_bgl_4.0_osx64.zip
%
% Mark Humphries, Mat Evans 28/2/2017

clear all; close all

fname = 'Lesmis.mat'; 
blnVizNet = 1;  % network visualisation - if MATLAB BGL installed, appropriate for platform:
fontsize = 6;

%% load data
load(['Results/Rejected_' fname],'Data','Rejection')

%% ready code for visualising networks
if blnVizNet
    % Traud Mucha Porter visualisation tools
    addpath('Traud_Mucha_Porter_CommunityVisualisation/');

    % needs MATLAB BGL Toolbox on your path - *change to your local path here*:
    if ismac
        bglpath = genpath('/Users/mqbssmhg/Dropbox/My Toolboxes/Graph_theory/matlab_bglOSX64/');  % generate path to local BGL and all its subdirectories
    else
        bglpath = genpath('C:\Users\mqbssmhg.DS\Dropbox\My Toolboxes\Graph_theory\matlab_bgl\');
    end
    % add to current MATLAB path
    addpath(bglpath); 
end


%% Image plot of A, with labels
k = sum(Data.A);
[~,iK] = sort(k,'descend');   % plot sorted into degree order
cmap = brewermap(10,'Greys');
figure; 
imagesc(Data.A(iK,iK)); colormap(cmap);
set(gca,'Ytick',1:length(Data.ixRetain));
set(gca,'Yticklabel',Data.nodelabels(iK,:),'Fontsize',fontsize);

%% compare eigenvalues of data and model
B = Data.A - Data.ExpA;  % modularity matrix using chosen null model
Edata = eig(B);
[fdata,xdata] = ecdf(Edata);
[fmodel,xmodel] = ecdf(Data.Emodel(:));
figure
stairs(xdata,fdata,'r'); hold on
stairs(xmodel,fmodel,'k')
xlabel('P')
ylabel('Eigenvalues')


%% visualise signal and noise parts
if blnVizNet
    n = size(Data.A,1);
    xynew = fruchterman_reingold_force_directed_layout(sparse(Data.A));
    syms = repmat('o',n,1);

    colors = repmat([0 0 0],n,1);
    colors(Rejection.ixNoise,:) = colors(Rejection.ixNoise,:) + 0.6;  % gray for noise 
    figure
    graphplot2D(xynew,Data.A,10,colors,syms,10);
    axis off
    allh = get(gca,'Children');
    set(allh,'MarkerEdgeColor',[0 0 0])
    set(allh,'LineWidth',0.2)
    title('Signal/Noise split')
end

% Add node labels to graph
for i = 1:length(Data.ixRetain); 
    txt(i) = text(xynew(i,1),xynew(i,2),Data.nodelabels(i,:),'Fontsize',fontsize);%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
end

%% Plot lowD projection of each node, with labels, sorted by magnitude
figure
[sorted_norms,SNIdx] = sort(Rejection.Difference.Norm); % SNIdx = Sorted Norm Index

stem(1:numel(Rejection.ixNoise),sorted_norms(1:numel(Rejection.ixNoise)));
hold all
stem(numel(Rejection.ixNoise)+1:numel(Rejection.Difference.Norm),sorted_norms(numel(Rejection.ixNoise)+1:end))
xlabel('Nodes')
ylabel('Projection (normalised to null model)')

% EDIT HERE: text labels 90 rotated: up for y>0; down for y<0
for i = 1:length(Data.A); 
    text(i,sorted_norms(i),Data.nodelabels(SNIdx(i),:),'Fontsize',fontsize);%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
end

%% visualise connected-signal network
if blnVizNet
    % visualise Aconnected
    colors = repmat([0 0 0],n,1)+0.6;
    colors(Data.ixConnectedSignal,:) = 0;  % black for signal, connected 
    figure
    graphplot2D(xynew,Data.A,10,colors,syms,10);
    axis off
    allh = get(gca,'Children');
    set(allh,'MarkerEdgeColor',[0 0 0])
    set(allh,'LineWidth',0.5)
    title('Connected Signal/Noise split')
end




