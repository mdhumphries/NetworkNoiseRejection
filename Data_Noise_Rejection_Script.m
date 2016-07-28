%% template script for applying complete work-flow to one data network, using one choice of null model
% data network: correlations between firing in Aplysia recording
% null model: weighted configuration model

clear all; close all

% analysis parameters
N = 20;        % repeats of permutation
alpha = 0.05;  % rejection region for noise
options.Weight = 'linear'; % 'linear' is default

% load data
load('Networks/Lesmis.mat');
A = full(Problem.A);

% get expected distribution of eigenvalues under null model (here, WCM)
% [Emodel,diagnostics] = expectedEigsUnd(A,N);
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
R = NodeRejection(B,Emodel(:),alpha,Vmodel,options); % N.B. also calls function to find projections

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

% [Cfull,Qmaxfull,Cconfull,Qcfull,Nfull,~] = allevsplitConTransitive(A);

%% plot sorted into group order
[srt,I] = sort(C,'ascend');
lines = [0; find(diff(srt)==1); numel(C)]+0.5;

% [srt,I] = sort(Cfull,'ascend');
% lines = [0; find(diff(srt)==1); numel(Cfull)]+0.5;

figure
% imagesc(Asignal(I,I);
imagesc(Aconnected(I,I));
% imagesc(A(I,I));

% draw outline box around each cluster
for i=2:numel(lines)
    line([lines(i-1) lines(i) lines(i) lines(i-1); lines(i) lines(i) lines(i-1) lines(i-1)],...
         [lines(i-1) lines(i-1) lines(i) lines(i); lines(i-1) lines(i) lines(i) lines(i-1)],...,
         'Color',[1 1 1],'LineWidth',1)
end
