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

if blnViz
    % Traud Mucha Porter visualisation tools
    addpath('Traud_Mucha_Porter_CommunityVisualisation/');

    % needs MATLAB BGL Toolbox on your path - change to your local path
    % here:
    bglpath = genpath('/Users/mqbssmhg/Dropbox/My Toolboxes/Graph_theory/matlab_bglOSX64/');  % generate path to local BGL and all its subdirectories
    
    % add to current MATLAB path
    addpath(bglpath); 
end

% analysis parameters
N = 100;        % repeats of permutation
alpha = 0.95;  % confidence interval on estimate of maxiumum eigenvalue for null model
options.Weight = 'linear'; % 'linear' is default
options.Weight = 'linear'; % 'linear' is default
options.Norm = 'L2'; % L2 is default

% load data
load('Networks/Lesmis.mat');
% load('Networks/dolphins.mat');
A = full(Problem.A);

% get expected distribution of eigenvalues under null model (here, WCM)
[Emodel,diagnostics,Vmodel] = WeightedConfigModel(A,N);

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

%% visualise signal and noise parts
if blnViz
    n = size(A,1);
    xynew = fruchterman_reingold_force_directed_layout(sparse(A));
    syms = repmat('o',n,1);

    colors = repmat([0 0 0],n,1);
    colors(R.ixNoise,:) = colors(R.ixNoise,:) + 0.6;  % gray for noise 
    figure
    graphplot2D(xynew,A,1,colors,syms,10);
    axis off
    allh = get(gca,'Children');
    set(allh,'MarkerEdgeColor',[0 0 0])
    set(allh,'LineWidth',0.2)
    title('Signal/Noise split')
end

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
title('Consensus clustering')

for i=1:numel(allC)
    CLou = allC{i}{1};  % Repeat#, Level of Hierarchy
    HL = plotClusterMap(Aconnected,CLou);
    title(['Louvain ' num2str(i)])
end