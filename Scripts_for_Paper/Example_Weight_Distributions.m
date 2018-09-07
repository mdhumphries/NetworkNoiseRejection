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
    stairs(Sparse(iN).xECDF,Sparse(iN).ECDF,'Color',[0.5 0.6 0.8],'Linewidth',0.4);

end
stairs(Data.xECDF,Data.ECDF,'k','Linewidth',2); hold on
set(gca,'Yscale','log')

% (2) as histogram
binW = [0 max(A(index))];

[Data.Hist,Data.Edges] = histcounts(A(index),'BinLimits',binW,'BinMethod','integers');
for iN = 1:pars.N
    Afull = poissonFull.A(:,:,iN);
    [Full(iN).Hist] = histcounts(Afull(index),'BinLimits',binW,'BinMethod','integers');
    Full(iN).Diff = Data.Hist - Full(iN).Hist;
    Asparse = poissonSparse.A(:,:,iN);
    [Sparse(iN).Hist] = histcounts(Asparse(index),'BinLimits',binW,'BinMethod','integers');
    Sparse(iN).Diff = Data.Hist - Sparse(iN).Hist;
end

centres = (Data.Edges(1:end-1) + Data.Edges(2:end))/2;
centres(1) = binW(1); centres(end) = binW(end);

% plot all histograms (skipping 0 entry)
figure
for iN = 1:pars.N
    stairs(centres(2:end),Full(iN).Hist(2:end),'Color',[0.8 0.6 0.5],'Linewidth',0.4); hold on
    stairs(centres(2:end),Sparse(iN).Hist(2:end),'Color',[0.5 0.6 0.8],'Linewidth',0.4); hold on
    
end
stairs(centres(2:end),Data.Hist(2:end),'k','Linewidth',2)

% set(gca,'Yscale','log','XLim',[-0.5 35],'YLim',[10^-1 10^4])

% plot difference in counts
figure
for iN = 1:pars.N
    stairs(centres,Full(iN).Diff,'Color',[0.8 0.6 0.5],'Linewidth',0.4); hold on
    stairs(centres,Sparse(iN).Diff,'Color',[0.5 0.6 0.8],'Linewidth',0.4); hold on 
end

% plot proportional difference in counts
figure
for iN = 1:pars.N
    stairs(centres,Full(iN).Diff ./ Data.Hist,'Color',[0.8 0.6 0.5],'Linewidth',0.4); hold on
    stairs(centres,Sparse(iN).Diff ./ Data.Hist,'Color',[0.5 0.6 0.8],'Linewidth',0.4); hold on 
end

save ../Results/WeightDistributionExample Data Full Sparse pars optionsModel binW
