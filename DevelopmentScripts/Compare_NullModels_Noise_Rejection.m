% template script for applying complete work-flow to one data network, using one choice of null model
% data network: correlations between firing in Aplysia recording
% null model: weighted configuration model
%
% Visualisations need:
% (1) Traud-Mucha-Porter toolbox (included in GitHub)
% (2) MATLAB BGL Toolbox:
%           Win32, Win64, Mac32, Linux: https://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/ 
%           Mac64: http://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/old/matlab_bgl_4.0_osx64.zip

clear all; 
close all
blnViz = 0;  % if MATLAB BGL installed, appropriate for platform:

if blnViz
    % Traud Mucha Porter visualisation tools
    addpath('Traud_Mucha_Porter_CommunityVisualisation/');

    % needs MATLAB BGL Toolbox on your path - change to your local path
    % here:
    % bglpath = genpath('/Users/mqbssmhg/Dropbox/My Toolboxes/Graph_theory/matlab_bglOSX64/');  % generate path to local BGL and all its subdirectories
    bglpath = genpath('C:\Users\lpzmdh\Dropbox\My Toolboxes\Graph_theory\matlab_bgl\');
    
    % add to current MATLAB path
    addpath(bglpath); 
end

addpath('../Network_Spectra_Functions/');

% analysis parameters
N = 100;        % repeats of permutation

%% LOAD DATA
% Example 1 Aplysia
% load Sep2711da01_40_120s_Sxy_Gaussian_1s
% A = Sxyall{1};
% Example 2
% load('lesmis.mat');
% A = full(Problem.A);
% Example 3
load('StarWarsNetworkEp1.mat');
% load('StarWarsNetworkEp2.mat');
% load('StarWarsNetworkEp3.mat');
A = full(StarWars.A);

%********************************
% get expected distribution of eigenvalues under null model (here, WCM)
% Choice of the models
tic
[Emodel{1},diagnostics{1},Vmodel{1}] = WeightedConfigModel(A,N,100);
rectime(1) = toc

tic
[Emodel{2},diagnostics{2},Vmodel{2}] = RndPoissonConfigModel(A,N);
rectime(2) = toc

tic
[Emodel{3},diagnostics{3},Vmodel{3}] = F_EmpiricalRndMat(A,N);
rectime(3) = toc
% decompose nodes into signal and noise
%*********************************
B = A - expectedA(A);  % modularity matrix using chosen null model

% compare data and model
Edata = eig(B);

color = ['r', 'b','k', 'g', 'm'];
% ###########################
% Quantify distribution of weights for: data, original configuration model
% & WCM
% WEIGHT DISTRIBUTION
figure()
for t = 1:length(Emodel);
    Aperm = [];
    for iN = 1:N;
        Aperm = [Aperm; nonzeros(diagnostics{t}(iN).Aperm)];
    end
    [f_perm{t},xi_perm{t}] = ksdensity(Aperm);

    plot(xi_perm{t},f_perm{t}./sum(f_perm{t})',color(t),'LineWidth',2); hold on
end
CM = expectedA(A);
[f_data,xi_data] = ksdensity(nonzeros(A));
[f_dataCM,xi_dataCM] = ksdensity(nonzeros(CM));
plot(xi_data,f_data./sum(f_data)',color(t+1),'LineWidth',2); hold on
plot(xi_dataCM,f_dataCM./sum(f_dataCM)',color(t+2),'LineWidth',2); hold on
legend('WCM','Poiss','Empir Rnd','Data','CM')
xlabel('Weight')
ylabel('Empirical PDF')
% ###########################
% Poisson and WCM asymptotically the same: over repeated generated null 
% models, show agreement for example network on degree distribution, 
% weight distribution, and distribution of maximum eigenvalue (latter can 
% only be shown over the multiple null models).
% DEGREE DISTRIBUTION
degreeEr = [];
degreeA = sum(A);
degree = [];
for t = 1:length(Emodel);
    deg = [];
    deger = [];
    for iN = 1:N;
        deg = [deg; sum(diagnostics{t}(iN).Aperm)'];
        deger = [deger; (sum(diagnostics{t}(iN).Aperm)-sum(A))'];
    end
    degree = [degree deg];
    degreeEr = [degreeEr deger];

end
figure()
subplot(121)
boxplot(degreeA,{'Data'}); hold on
ylabel('Degree distribution')
title('Original network')
ylim([min(min(degree)) max(max(degree))])
subplot(122)
boxplot(degree,{'WCM','Poiss','Empir rnd'}); hold on
title('Null model')
ylim([min(min(degree)) max(max(degree))])
figure()
boxplot(degreeEr,{'WCM','Poiss','Empir rnd'}); hold on
ylabel('Delta Degree distribution')
% MAX EIGENVALUE DISTRIBUTION
for t = 1:length(Emodel);
    MaxEig = max(Emodel{t});
    % test for normal distribution of the maximum eigenvalues distribution.
    % If h=0 we cannot reject the nul hypothesis that the MaxEig is
    % normally distributed.
    [h, p] = chi2gof(MaxEig);
%    % if uncomment the following lines, estimate the normal parameters only
%    if the condition of normality is satisfied.
%     if h==0;
        [pd] = fitdist(MaxEig','Normal');
%     else 
%         dist('Maximum Eigenvalues distribution is not normal')
%         return
%     end

    % some visualization of CDF for all the max Eigenv of the null model and
    % from a random gaussian distribution with mean and std estimated from the
    % null model
    [feigdata,xeigdata] = ecdf(MaxEig);
    [fGaussmodel, xGaussmodel] = ecdf(pd.mu+pd.sigma*randn(1000000,1));

    figure(30)
    stairs(xeigdata,feigdata,'r'); hold on
    stairs(xGaussmodel,fGaussmodel,'k'); hold on
    legend('Data','GaussModel')
    title('comparison of null models Max eigenvalues distribution')
    xlabel('Max Eigenvalue')
    ylabel('CDF')
    
    % Monte Carlo simulation. Generate 1 million sample of size N (#
    % permutation) and determine the maximum for each (total 1 million eig).
    % Then randomly sample from data, with replacement
    npoints = 1000000; 
    sample = pd.mu+pd.sigma*randn(length(MaxEig),npoints);
    SAMPLE = sample(:);

    %Compare the distribution of the maximum eigenvalue between Emodel and the
    %Monte Carlo distribution
    [binc,binrange] = ksdensity(MaxEig);
    [bincMC,binrangeS] = ksdensity(SAMPLE);
    [f_max,x_maxeig] = ksdensity(MaxEig);
    ind{t} = find(Edata>max(MaxEig));
    
    figure(40)
    plot(x_maxeig,f_max./sum(f_max)',color(t),'LineWidth',2); hold on
%     plot(Edata(ind{t}),zeros(length(ind{t}),1),color(t),'Marker','o','MarkerFaceColor','r'); hold on
    ylabel('PDF')
    xlabel('Max Eigenvalues')
end
figure(40)
legend('WCM','Poiss','Empir rnd')

% EIGENVALUE DISTRIBUTION
figure()
for t = 1:length(Emodel);
    [f_eig,x_eig] = ksdensity(Emodel{t}(:));
    plot(x_eig,f_eig./sum(f_eig)',color(t),'LineWidth',2); hold on
end
legend('WCM','Poiss','Empir Rnd')
ylabel('PDF')
xlabel('Eigenvalues')
    
