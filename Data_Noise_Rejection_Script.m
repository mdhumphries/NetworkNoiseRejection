%% template script for applying complete work-flow to one data network, using one choice of null model
% data network: correlations between firing in Aplysia recording
% null model: weighted configuration model

clear all; close all

% analysis parameters
N = 100;        % repeats of permutation
alpha = 0.05;  % rejection region for noise
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
[Dspace,Dn] = LowDSpace(B,Emodel(:),alpha); % to just obtain low-dimensional projection

% node rejection within low-dimensional projection
R = NodeRejection(B,Emodel,alpha,Vmodel,options); % N.B. also calls function to find projections

% new signal matrix
Asignal = A(R.ixSignal,R.ixSignal);

%% analyse new signal matrix

% % first: remove any nodes without connections
kAsignal = sum(Asignal>0);
ixConnectedSignal = R.ixSignal(kAsignal > 1);  % more than 1 link
Aconnected = A(ixConnectedSignal,ixConnectedSignal);  % subset of original matrix

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