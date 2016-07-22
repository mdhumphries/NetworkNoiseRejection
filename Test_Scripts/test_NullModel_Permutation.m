%% script to illustrate the permutation testing script

clear all; close all

addpath('../')

load Sep2711da01_40_120s_Sxy_Gaussian_1s

A = Sxyall{1};
N = 20;  % repeats of permutation

% do permutations
[allV,diagnostics] = expectedEigsUnd(A,N);

%% get bounds on configuration model eigenvalues
bnds95 = prctile(allV,[2.5,97.5]); % 95% confidence interval on eigenvalue distribution for null model
bnds99 = prctile(allV,[0.5,99.5]);
[h,x] = hist(allV,20);

figure
bar(x,h,1,'Facecolor',[1 1 1])
m = max(get(gca,'YLim'));
line([bnds95(1) bnds95(1)],[0 m],'Color','b')
line([bnds95(2) bnds95(2)],[0 m],'Color','b')
line([bnds99(1) bnds99(1)],[0 m],'Color','r')
line([bnds99(2) bnds99(2)],[0 m],'Color','r')
xlabel('Eigenvalue')
ylabel('Frequency')

%% compare to data eigenvalues

[V,D] = eig(A);
egs = sort(diag(D),'descend');

% number of communities
Ncommunities95 = find(egs >= bnds95(2));
Ncommunities99 = find(egs >= bnds99(2));

% number of "noise" nodes
Nnoise95 = find(egs > bnds95(1) & egs < bnds95(2));
Nnoise99 = find(egs > bnds99(1) & egs < bnds99(2));





