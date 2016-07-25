%% script to illustrate the permutation testing script

clear all; close all

addpath('../')
N = 20;  % repeats of permutation


%% no communities example
A = rand(100)./2;
A = A + A';  % undirected, weighted network in [0,1]
A(eye(100)==1) = 0;
% A(A < 0.1) = 0; 

% do permutations
[allV,diagnostics] = expectedEigsUnd(A,N);

% get bounds on configuration model eigenvalues
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
title('Random weighted graph')

Bdata = A - expectedA(A);
[V,D] = eig(Bdata);
egs = sort(diag(D),'descend');

% number of communities
Ncommunities95 = find(egs >= bnds95(2));
Ncommunities99 = find(egs >= bnds99(2));

% number of "noise" nodes
Nnoise95 = find(egs > bnds95(1) & egs < bnds95(2));
Nnoise99 = find(egs > bnds99(1) & egs < bnds99(2));

% compare eigenvalues
[f,x] = ecdf(egs);
[fperm,xperm] = ecdf(allV);
figure
stairs(x,f,'r'); hold on
stairs(xperm,fperm,'k')

%% configuration model example

%% data example
load Sep2711da01_40_120s_Sxy_Gaussian_1s

A = Sxyall{1};

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
title('Data weighted graph')

% compare to M-P bounds
B = A - expectedA(A);
s = std(B(:));
c = 1; % ratio of MxN matrix
MP_low =(s^2)*(1-sqrt(c))^2;
MP_high =(s^2)*(1+sqrt(c))^2;

line([MP_low MP_low],[0 m],'Color',[0.6 0.6 0.6])
line([MP_high MP_high],[0 m],'Color',[0.6 0.6 0.6])


%% compare to data eigenvalues

[V,D] = eig(A);
egs = sort(diag(D),'descend');

% number of communities
Ncommunities95 = find(egs >= bnds95(2));
Ncommunities99 = find(egs >= bnds99(2));

% number of "noise" nodes
Nnoise95 = find(egs > bnds95(1) & egs < bnds95(2));
Nnoise99 = find(egs > bnds99(1) & egs < bnds99(2));





