%% template script for applying complete work-flow to one data network, using one choice of null model
% data network: correlations between firing in Aplysia recording
% null model: weighted configuration model

clear all; close all

% analysis parameters
N = 20;        % repeats of permutation
alpha = 0.05;  % rejection region for noise

% load data
load('Test_Scripts/Sep2711da01_40_120s_Sxy_Gaussian_1s')
A = Sxyall{1};

% get expected distribution of eigenvalues under null model (here, WCM)
[allV,diagnostics] = expectedEigsUnd(A,N);

% decompose nodes into signal and noise
B = A - expectedA(A);  % modularity matrix using null model
D = NodeRejection(B,allV,alpha);

% new signal matrix
Asignal = A(D.ixSignal,D.ixSignal);

%% analyse new signal matrix

% consensus modularity
[C,Qmax,Ccon,Qc,N,Q] = allevsplitConTransitive(Asignal);

%% plot sorted into group order
[srt,I] = sort(Ccon,'ascend');
lines = [0; find(diff(srt)==1); numel(C)]+0.5;
figure
imagesc(Asignal(I,I))
% draw outline box around each cluster
for i=2:numel(lines)
    line([lines(i-1) lines(i) lines(i) lines(i-1); lines(i) lines(i) lines(i-1) lines(i-1)],...
         [lines(i-1) lines(i-1) lines(i) lines(i); lines(i-1) lines(i) lines(i) lines(i-1)],...,
         'Color',[1 1 1],'LineWidth',1)
end
