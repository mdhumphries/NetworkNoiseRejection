%% template script for applying complete work-flow to one data network, using one choice of null model
% data network: correlations between firing in Aplysia recording
% null model: weighted configuration model
% Mark Humphries 28/2/2017

clear all; close all

% network to analyse
fname = 'LesMis.mat'; 

% analysis parameters
pars.N = 100;           % repeats of permutation
pars.alpha = 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
pars.Model = 'Poiss';   % or 'WCM' . % which null model
pars.C = 1;             % conversion factor for real-valued weights (set=1 for integers)

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble? 
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?

% NodeRejection options
optionsReject.Weight = 'linear'; % 'linear' is default
optionsReject.Norm = 'L2';       % L2 is default

%% load data-file
load(['Networks/' fname]); 

if strfind(fname,'StarWars')
    A = StarWars.A;
    nodelabels = StarWars.Nodes;
    Data.nodelabels = nodelabels';
elseif strfind(fname,'cosyne')
    A = adjMatrix;
    m = cellfun('length',cosyneData.authorHash);
    nodelabels = [];
    for i = 1:numel(cosyneData.authorHash)
        l = numel(cosyneData.authorHash{i});
        nodelabels = [nodelabels; cosyneData.authorHash{i} blanks(max(m) - l)];
    end
    Data.nodelabels = nodelabels;
else
    A = full(Problem.A);
    % Generate node labels for later visualisation to work
    Data.nodelabels = Problem.aux.nodename;
end

% clean-up A, and store as basis for all further analysis
[Data.A,Data.ixRetain] = prep_A(A);

% % SBM generation
% A_SBM = test_noise_rejection_planted_noise(50,2,'low',0.2);
% A = A_SBM.adjacency;
% nodelabels = num2str(A_SBM.membership);
% A = round(A);

%% get expected distribution of eigenvalues under null model (here, WCM)

switch pars.Model
    case 'Poiss'
        [Data.Emodel,diagnostics,Vmodel,Data.ExpA] = RndPoissonConfigModel(Data.A,pars.N,pars.C,optionsModel);
    case 'WCM'
        [Data.Emodel,diagnostics,Vmodel,Data.ExpA] = WeightedConfigModel(Data.A,pars.N,pars.C,optionsModel);
    otherwise
        error('Unrecognised null model specified')
end

%% decompose nodes into signal and noise
B = Data.A - Data.ExpA;  % modularity matrix using chosen null model

% find low-dimensional projection
[Data.Dspace,~,Data.Dn,Data.EigEst] = LowDSpace(B,Data.Emodel,pars.alpha); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors

% node rejection within low-dimensional projection
Rejection = NodeRejection(B,Data.Emodel,pars.alpha,Vmodel,optionsReject); % N.B. also calls LowDSpace function to find projections

% new signal matrix
Data.Asignal = Data.A(Rejection.ixSignal,Rejection.ixSignal);

% connected signal matrix
kAsignal = sum(Data.Asignal>0);
Data.ixConnectedSignal = Rejection.ixSignal(kAsignal > 1);  % more than 1 link
Data.Aconnected = Data.A(Data.ixConnectedSignal,Data.ixConnectedSignal);  % subset of original matrix

%% compare to standard configuration model
[Control.Emodel,diagnostics,Vmodel] = RndPoissonConfigModel(Data.A,pars.N,pars.C);
Control.P = expectedA(Data.A);

B = Data.A - Control.P;
[Control.Dspace,~,Control.Dn,Control.EigEst] = LowDSpace(B,Control.Emodel,pars.alpha); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors




%% save
save(['Results/Rejected_' fname],'Rejection','Data','Control','pars','optionsModel','optionsReject')
