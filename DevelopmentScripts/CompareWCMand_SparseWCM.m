% script to test all null models...
clearvars

addpath ../Helper_Functions/
addpath ../Network_Spectra_Functions/
addpath ../Network_Analysis_Functions/

load '../Networks/LesMis';
A = full(Problem.A);
pars.C = 1;             % conversion factor for real-valued weights (set=1 for integers)

%load Sep2711da01_40_120s_Sxy_Gaussian_1s.mat
% A = Sxyall{1};
% pars.C = 100;

 % analysis parameters
pars.N = 100;           % repeats of sampling
% pars.alpha = 0; %0.95; % 0.95; % 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
pars.I = 0;      % interval: set to 0 for mean
pars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble? 
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?


%% run each null model...

[linkFull.E,linkFull.D,linkFull.V,linkFull.A] = linkFullWCM(A,pars.N,pars.C);

[poissonFull.E,poissonFull.D,poissonFull.V,poissonFull.A] = poissonFullWCM(A,pars.N,pars.C);

[linkSparse.E,linkSparse.D,linkSparse.V,linksparse.ExpA,linkSparse.A] = linkSparseWCM(A,pars.N,pars.C,optionsModel);

[poissonSparse.E,poissonSparse.D,poissonSparse.V,poissonsparse.ExpA,poissonSparse.A] = poissonSparseWCM(A,pars.N,pars.C,optionsModel);


%% check output
plotY = cell(4,1);
plotY{1} = [linkFull.D(:).dStotal];
plotY{2} = [poissonFull.D(:).dStotal];
plotY{3} = [linkSparse.D(:).dStotal];
plotY{4} = [poissonSparse.D(:).dStotal];

figure
UnpairedUnivariateScatterPlots(gca,plotY,'strX',{'linkF','PoissF','linkS','PoissS'});
ylabel('Difference in total strength')

plotY{1} = [linkFull.D(:).dDensity];
plotY{2} = [poissonFull.D(:).dDensity];
plotY{3} = [linkSparse.D(:).dDensity];
plotY{4} = [poissonSparse.D(:).dDensity];

figure
UnpairedUnivariateScatterPlots(gca,plotY,'strX',{'linkF','PoissF','linkS','PoissS'});
ylabel('Difference in network density')

