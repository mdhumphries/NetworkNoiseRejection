%% template script for applying complete work-flow to one data network, using one choice of null model
% data network: correlations between firing in Aplysia recording
% null model: weighted configuration model
% Mark Humphries 28/2/2017

clear all; close all

% network to analyse
% fname = 'Allen_Gene_Leaf'; 
fname = 'LesMis'; 


% analysis parameters
pars.N = 100;           % repeats of permutation
% pars.alpha = 0; %0.95; % 0.95; % 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
pars.I = 0;      % interval 
pars.Model = 'Poiss';   % or 'WCM' . % which null model
pars.C = 1;             % conversion factor for real-valued weights (set=1 for integers)
pars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble? 
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?

% NodeRejection options
optionsReject.Weight = 'linear'; % 'linear' is default
optionsReject.Norm = 'L2';       % L2 is default
optionsReject.Interval = 'CI';

%% load data-file
load(['Networks/' fname]); 

if strfind(fname,'StarWars')
    A = StarWars.A;
    nodelabels = StarWars.Nodes;
    nodelabels = nodelabels';
elseif strfind(fname,'cosyne')
    A = adjMatrix;
    m = cellfun('length',cosyneData.authorHash);
    nodelabels = [];
    for i = 1:numel(cosyneData.authorHash)
        l = numel(cosyneData.authorHash{i});
        nodelabels = [nodelabels; cosyneData.authorHash{i} blanks(max(m) - l)];
    end
    nodelabels = nodelabels;
elseif strfind(fname,'CosyneYear')
    A = adjMatrix;
    nodelabels = nodelabel;
elseif exist('Problem')
    A = full(Problem.A);
    % Generate node labels for later visualisation to work
    nodelabels = Problem.aux.nodename;
end

% make undirected if necessary
A = (A + A') / 2; % make undirected

% clean-up A, get largest component, and store as basis for all further analysis
% all indices are with reference to Data.A
[Data.A,Data.ixRetain,Data.Comps,Data.CompSizes] = prep_A(A);
Data.nodelabels = nodelabels(Data.ixRetain,:);

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
[Data.Dspace,~,Data.Dn,Data.EigEst,Data.Nspace,~,Data.Dneg,Data.NEigEst] = LowDSpace(B,Data.Emodel,pars.I); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors

% compute dimensions based on just positive eigenvalues
egs = eig(B);  % eigenspectra of data modularity matrix
egs = sort(egs,'descend'); % sort eigenvalues into descending order 
Data.PosDn = sum(egs > pars.eg_min);

% node rejection within low-dimensional projection
Rejection = NodeRejection(B,Data.Emodel,pars.I,Vmodel,optionsReject); % N.B. also calls LowDSpace function to find projections

% new signal matrix
Data.Asignal = Data.A(Rejection.ixSignal,Rejection.ixSignal);

% connected signal matrix: find largest component, and use that - store
% others
[Data.Asignal_comp,ixRetain,Data.SignalComps,Data.SignalComp_sizes] = prep_A(Data.Asignal); 
Data.ixSignal_comp = Rejection.ixSignal(ixRetain);  % original node indices

% and then strip out leaves - nodes with single links
K = sum(Data.Asignal_comp);
ixLeaves = find(K==1); ixKeep = find(K > 1);

Data.ixSignal_Final = Data.ixSignal_comp(ixKeep);
Data.ixSignal_Leaves = Data.ixSignal_comp(ixLeaves);
Data.Asignal_final = Data.Asignal_comp(ixKeep,ixKeep);

%% compare to standard configuration model
[Control.Emodel,diagnostics,Vmodel] = RndPoissonConfigModel(Data.A,pars.N,pars.C);
Control.P = expectedA(Data.A);

B = Data.A - Control.P;

% compute groups based on just positive eigenvalues
egs = eig(B);  % eigenspectra of data modularity matrix
egs = sort(egs,'descend'); % sort eigenvalues into descending order 
Control.PosDn = sum(egs > pars.eg_min);

% compute groups based on estimated bounds
[Control.Dspace,~,Control.Dn,Control.EigEst] = LowDSpace(B,Control.Emodel,pars.I); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors




%% save
save(['Results/Rejected_' fname],'Rejection','Data','Control','pars','optionsModel','optionsReject')
