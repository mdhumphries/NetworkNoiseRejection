%% script to assess performance of group detection as tuning network from random to modular 
% Compare:
%       Communities detected by community detection algorithms (using Q >
%       0)
%       VS
%       Communities detected by spectral rejection
%
% Here in script 01: create networks and detect with spectral rejection
%
% Change log:
% 5/6/2018: initial version
%
% Mark Humphries

clear all; close all;
addpath('../SyntheticModel/');
addpath('../Network_Spectra_Functions/');

fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';

%% fixed parameters of synthetic model
% Model.N = [200,75,25];  % size of modules
Model.N = [100,100,100];  % size of modules

Model.P.between = 0.05;   
Model.alpha = 0;    % weights follow links exactly

% stength distribution parameters
Model.Spar.distribution =  'Poisson';  % type of distribution
Model.Spar.a = 50;                    % scale: in addition to existing edges
Model.Spar.b = 1;                     % spread

nBatch = 100;

%% range of parameters
Model.P_of_within = [0.05:0.025:0.2];

%% rejection and clustering parameters
rejectionpars.N = 100;           % repeats of permutation
rejectionpars.alpha = 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
rejectionpars.Model = 'Poiss';   % or 'WCM' . % which null model
rejectionpars.C = 1;            % conversion factor for real-valued weights (set=1 for integers)
rejectionpars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble? 
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?

% NodeRejection options
optionsReject.Weight = 'linear'; % 'linear' is default
optionsReject.Norm = 'L2';       % L2 is default


%% loop: make models and do rejection etc
n = sum(Model.N);
P.between = Model.P.between;

% group membership
Nsum = [0 cumsum(Model.N)];
group_membership = zeros(n,1);
for iG = 1:numel(Nsum)-1
    group_membership(Nsum(iG)+1:Nsum(iG+1)) = iG;
end 


for iP = 1:numel(Model.P_of_within)
    % assign parameters   
    P.in = Model.P_of_within(iP);

    for iB = 1:nBatch
        disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_within)) '; batch: ' num2str(iB)])
        tic
        %% make model
        A = wire_edges(Model.N,P);  % make adjacency matrix
        S = sample_strength(n,Model.Spar); % sample strength distribution according to current rules            
        W = weight_edges(A,Model.N,S,Model.alpha,P); % use Poisson generative model to create Weight matrix

        %% do spectral rejection on model
        [Data,Rejection,Control] = reject_the_noise(W,rejectionpars,optionsModel,optionsReject);        

        % store relevant network information for clustering
        Network(iP,iB).W = W;
        Network(iP,iB).ExpW = Data.ExpA;

        % RESULTS:
        % number of groups recovered by spectra vs same count from other
        % approaches

        Results.SpectraWCMGroups(iP,iB) = Data.Dn+1;
        Results.SpectraConfigGroups(iP,iB) = Control.Dn+1;
        Results.PosEigWCMGroups(iP,iB) = Data.PosDn+1;
        Results.PosEigConfigGroups(iP,iB) = Control.PosDn+1;

        Results.Time(iP,iB) = toc;
    end
end

%% process results
Results.ProportionModular.SpectraWCM = sum(Results.SpectraWCMGroups > 1,2);
Results.ProportionModular.SpectraConfig = sum(Results.SpectraConfigGroups > 1,2);
Results.ProportionModular.PosEigWCM = sum(Results.PosEigWCMGroups > 1,2);
Results.ProportionModular.PosEigConfig = sum(Results.PosEigConfigGroups > 1,2);

    
%% save
fname = datestr(now,30);

if all(Model.N == Model.N(1))
    strName = 'Equal';
else
    strName = 'Unequal';
end

save([fpath 'P_rejection_Synthetic' strName '_NoNoise_' fname],'Results','Network','Model','group_membership','rejectionpars','optionsModel','optionsReject')


