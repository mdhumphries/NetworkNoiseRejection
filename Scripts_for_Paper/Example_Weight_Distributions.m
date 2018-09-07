%% example weight distributions

clearvars; close all

addpath ../Helper_Functions/
addpath ../Network_Spectra_Functions/
addpath ../Network_Analysis_Functions/

% load '../Networks/celegansneural.mat';
% A = full(Problem.A);
% A = (A + A') /2;
% pars.C = 2;             % conversion factor for real-valued weights (set=1 for integers)

load '../Networks/lesmis.mat';
A = full(Problem.A);
pars.C = 1;             % conversion factor for real-valued weights (set=1 for integers)

 % analysis parameters
pars.N = 100;           % repeats of sampling
% pars.alpha = 0; %0.95; % 0.95; % 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
pars.I = 0;      % interval: set to 0 for mean
pars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble? 
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?


%% run each null model...

% [linkFull.E,linkFull.D,linkFull.V,linkFull.A] = linkFullWCM(A,pars.N,pars.C);

[poissonFull.E,poissonFull.D,poissonFull.V,poissonFull.A] = poissonFullWCM(A,pars.N,pars.C);

% [linkSparse.E,linkSparse.D,linkSparse.V,linksparse.ExpA,linkSparse.A] = linkSparseWCM(A,pars.N,pars.C,optionsModel);

[poissonSparse.E,poissonSparse.D,poissonSparse.V,poissonsparse.ExpA,poissonSparse.A] = poissonSparseWCM(A,pars.N,pars.C,optionsModel);


%% weight distributions

% just use unique weights, and not diagonal
index = find(tril(ones(size(A)), -1));

% (1) as CDF
[Data.ECDF,Data.xECDF] = ecdf(A(index));

for iN = 1:pars.N
    Afull = poissonFull.A(:,:,iN);
    [Full(iN).ECDF,Full(iN).xECDF] = ecdf(Afull(index));
    Asparse = poissonSparse.A(:,:,iN);
    [Sparse(iN).ECDF,Sparse(iN).xECDF] = ecdf(Asparse(index));
end

figure
for iN = 1:pars.N
    stairs(Full(iN).xECDF,Full(iN).ECDF,'Color',[0.8 0.6 0.5],'Linewidth',0.4); hold on
    % stairs(Sparse(iN).xECDF,Sparse(iN).ECDF,'Color',[0.5 0.6 0.8],'Linewidth',0.4);

end
stairs(Data.xECDF,Data.ECDF,'k'); hold on
set(gca,'Yscale','log')