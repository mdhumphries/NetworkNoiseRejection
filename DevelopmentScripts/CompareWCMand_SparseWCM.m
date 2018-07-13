% script to test all null models...
 clearvars
 
 load '../Networks/LesMis';
 
 % analysis parameters
pars.N = 100;           % repeats of sampling
% pars.alpha = 0; %0.95; % 0.95; % 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
pars.I = 0;      % interval: set to 0 for mean
pars.C = 1;             % conversion factor for real-valued weights (set=1 for integers)
pars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble? 
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?

%% stats of original network

s = sum(Problem.A);
S = sum(s)/2;

density = sum(Problem.A > 2);

%% run each null model...

[linkFull.E,linkFull.D,linkFull.V,linkFull.A] = linkFullWCM(Problem.A,pars.N);

[linkSparse.E,linkSparse.D,linkSparse.V,linksparse.ExpA,linkSparse.A] = linkSparseWCM(Problem.A,pars.N,[],optionsModel);



