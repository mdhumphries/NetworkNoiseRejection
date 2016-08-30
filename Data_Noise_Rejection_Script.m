%% template script for applying complete work-flow to one data network, using one choice of null model
% data network: correlations between firing in Aplysia recording
% null model: weighted configuration model
%
% Visualisations need:
% (1) Traud-Mucha-Porter toolbox (included in GitHub)
% (2) MATLAB BGL Toolbox:
%           Win32, Win64, Mac32, Linux: https://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/ 
%           Mac64: http://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/old/matlab_bgl_4.0_osx64.zip

clear all; close all
blnViz = 1;  % if MATLAB BGL installed, appropriate for platform:

% if blnViz
%     % Traud Mucha Porter visualisation tools
%     addpath('Traud_Mucha_Porter_CommunityVisualisation/');
% 
%     % needs MATLAB BGL Toolbox on your path - change to your local path
%     % here:
%     % bglpath = genpath('/Users/mqbssmhg/Dropbox/My Toolboxes/Graph_theory/matlab_bglOSX64/');  % generate path to local BGL and all its subdirectories
%     bglpath = genpath('C:\Users\mqbssmhg.DS\Dropbox\My Toolboxes\Graph_theory\matlab_bgl\');
%     
%     % add to current MATLAB path
%     addpath(bglpath); 
% end


% analysis parameters
N = 100;        % repeats of permutation
alpha = 0.95;  % confidence interval on estimate of maxiumum eigenvalue for null model
options.Weight = 'linear'; % 'linear' is default
options.Norm = 'L2'; % L2 is default

% % load data
% load('Networks/Lesmis.mat');
% % load('Networks/dolphins.mat');
% A = full(Problem.A);
% % Generate node labels for later visualisation to work
% nodelabels = Problem.aux.nodename;

% SBM generation
A_SBM = test_noise_rejection_planted_noise(50,2,'low',0.2);
A = A_SBM.adjacency;
nodelabels = num2str(A_SBM.membership);


tic
A = round(A);
% get expected distribution of eigenvalues under null model (here, WCM)
[Emodel,diagnostics,Vmodel] = WeightedConfigModel(A,N,1);


% decompose nodes into signal and noise
B = A - expectedA(A);  % modularity matrix using chosen null model

% compare data and model
Edata = eig(B);
[fdata,xdata] = ecdf(Edata);
[fmodel,xmodel] = ecdf(Emodel(:));
figure
stairs(xdata,fdata,'r'); hold on
stairs(xmodel,fmodel,'k')
xlabel('P')
ylabel('Eigenvalues')

% find low-dimensional projection
[Dspace,Ix,Dn,EigEst] = LowDSpace(B,Emodel,alpha); % to just obtain low-dimensional projection

% node rejection within low-dimensional projection
R = NodeRejection(B,Emodel,alpha,Vmodel,options); % N.B. also calls function to find projections

% new signal matrix
Asignal = A(R.ixSignal,R.ixSignal);
toc

%% visualise signal and noise parts
if blnViz
    n = size(A,1);
    xynew = fruchterman_reingold_force_directed_layout(sparse(A));
    syms = repmat('o',n,1);

    colors = repmat([0 0 0],n,1);
    colors(R.ixNoise,:) = colors(R.ixNoise,:) + 0.6;  % gray for noise 
    figure
    graphplot2D(xynew,A,10,colors,syms,10);
    axis off
    allh = get(gca,'Children');
    set(allh,'MarkerEdgeColor',[0 0 0])
    set(allh,'LineWidth',0.2)
    title('Signal/Noise split')
end

%% Add node labels to graph
for i = 1:length(A); 
    txt(i) = text(xynew(i,1),xynew(i,2),nodelabels(i,:));%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
end

%% Remove node labels
for i = 1:length(A); txt(i).Color = 'none';end

%% Plot lowD projection of each node, with labels, sorted by magnitude
figure
[sorted_norms,SNIdx] = sort(R.Difference.Norm); % SNIdx = Sorted Norm Index

stem(1:numel(R.ixNoise),sorted_norms(1:numel(R.ixNoise)));
hold all
stem(numel(R.ixNoise)+1:numel(R.Difference.Norm),sorted_norms(numel(R.ixNoise)+1:end))

for i = 1:length(A); 
    text(i,sorted_norms(i),nodelabels(SNIdx(i),:));%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
end

% Unsorted version
% figure
% stem(R.ixSignal,R.Difference.Norm(R.ixSignal))
% hold all
% stem(R.ixNoise,R.Difference.Norm(R.ixNoise))
% 
% for i = 1:length(A); 
%     text(i,R.Difference.Norm(i),Problem.aux.nodename(i,:));%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
% end

%% analyse new signal matrix

% % first: remove any nodes without connections
kAsignal = sum(Asignal>0);
ixConnectedSignal = R.ixSignal(kAsignal > 1);  % more than 1 link
Aconnected = A(ixConnectedSignal,ixConnectedSignal);  % subset of original matrix

if blnViz
    % visualise Aconnected
    colors = repmat([0 0 0],n,1)+0.6;
    colors(ixConnectedSignal,:) = 0;  % black for signal, connected 
    figure
    graphplot2D(xynew,A,1,colors,syms,10);
    axis off
    allh = get(gca,'Children');
    set(allh,'MarkerEdgeColor',[0 0 0])
    set(allh,'LineWidth',0.5)
    title('Connected Signal/Noise split')
end

% consensus modularity
% [C,Qmax,Ccon,Qc,N,Q] = allevsplitConTransitive(Asignal);
[C,Qmax,Ccon,Qc,N,~] = allevsplitConTransitive(Aconnected);
% [C2,Qmax2,Ccon2,Qc2,N,~] = allevsplitConTransitive(Aconnected);

% [Cfull,Qmaxfull,Cconfull,Qcfull,Nfull,~] = allevsplitConTransitive(A);

% Louvain algorithm
[allC,allQ,allCn,allIters] = LouvainCommunityUDnondeterm(Aconnected,5,1);  % run 5 times; return 1st level of hierarchy only

%% plot sorted into group order

H = plotClusterMap(Aconnected,Ccon);
title('Consensus clustering');
% Add node labels
numConnected = length(Aconnected);
[srt,I] = sort(Ccon,'ascend');
set(gca,'Xtick',1:numConnected);
set(gca,'Xticklabel',nodelabels(I,:));
set(gca,'XTickLabelRotation',90);


for i=1:numel(allC)
    CLou = allC{i}{1};  % Repeat#, Level of Hierarchy
    HL = plotClusterMap(Aconnected,CLou);
    title(['Louvain ' num2str(i)]);
    % Add node labels
    [srt,I] = sort(CLou,'ascend');
    set(gca,'Xtick',1:numConnected);
    set(gca,'Xticklabel',nodelabels(I,:));
    set(gca,'XTickLabelRotation',90);
end