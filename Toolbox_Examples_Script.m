%% template script for applying complete work-flow to one data network
%
% null model: weighted configuration model
%
% 06/08/2018 : first version
% Mark Humphries 

clearvars; close all

addpath('Network_Spectra_Functions/')
addpath('Network_Analysis_Functions/')
addpath('Helper_Functions/')

% network to analyse
fname = 'LesMis'; 

% analysis parameters: weight conversion is set dynamically, see below
pars.N = 100;           % repeats of permutation
% pars.alpha = 0; %0.95; % 0.95; % 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
pars.I = 0;      % interval: set to 0 for mean
pars.Model = 'Poiss';   % or 'Link'; % which version of sparse WCM null model
pars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble? 
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?

% NodeRejection options
optionsReject.Weight = 'linear'; % 'linear' is default
optionsReject.Norm = 'L2';       % L2 is default
optionsReject.Interval = 'CI';   % rejection if fall within confidence interval - but only if pars.I is not 0!

%% load data-file
load(['Networks/' fname]); 

% deal with individual differences in data formats:
% A : full weight matrix
% nodelabels : a string array, one string per row
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
    if isfield(Problem,'aux')
        nodelabels = Problem.aux.nodename;
    else
        nodelabels = string(1:size(A,1))';
    end
end

% make undirected if necessary
A = (A + A') / 2; % make undirected

% clean-up A, get largest component, and store as basis for all further analysis
% all indices are with reference to Data.A
[Data.A,Data.ixRetain,Data.Comps,Data.CompSizes] = prep_A(A);
Data.nodelabels = nodelabels(Data.ixRetain,:);   % update the node labels

%% set conversion parameter
if all(Data.A(Data.A>0) == 1) || ~any(rem(Data.A(:),1)) % binary or integers
    pars.C = 1;
else % has real values
    switch fname
        case {'celegansneural','polblogs'}
            pars.C = 2;  % has just full and half counts (1,1.5, etc)
        otherwise
            % assume is continuous real-valued
            pars.C = 100;
    end
end


%% get expected distribution of eigenvalues under null model (here, sparse WCM)
% functions are also available for computing these using the full weighted
% configuration model
switch pars.Model
    case 'Poiss'
        [Data.Emodel,diagnostics,Vmodel,Data.ExpA] = poissonSparseWCM(Data.A,pars.N,pars.C,optionsModel);
    case 'Link'
        [Data.Emodel,diagnostics,Vmodel,Data.ExpA] = linkSparseWCM(Data.A,pars.N,pars.C,optionsModel);
    otherwise
        error('Unrecognised null model specified')
end

keyboard

%% decompose nodes into signal and noise
B = Data.A - Data.ExpA;  % modularity matrix using chosen null model

% find low-dimensional projection: this is an optional step
% we call this here to demonstrate its useful output; it is called internally by NodeRejection
[Data.Dspace,~,Data.Dn,Data.EigEst,Data.Nspace,~,Data.Dneg,Data.NEigEst] = LowDSpace(B,Data.Emodel,pars.I); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors

% node rejection within low-dimensional projection
Rejection = NodeRejection(B,Data.Emodel,pars.I,Vmodel,optionsReject); % N.B. also calls LowDSpace function to find projections

% new signal matrix
Data.Asignal = Data.A(Rejection.ixSignal,Rejection.ixSignal);

% connected signal matrix: find largest component, and use that
[Data.Asignal_comp,ixRetain,Data.SignalComps,Data.SignalComp_sizes] = prep_A(Data.Asignal); 
Data.ixSignal_comp = Rejection.ixSignal(ixRetain);  % original node indices

% and then strip out leaves - nodes with single links
K = sum(Data.Asignal_comp);
ixLeaves = find(K==1); ixKeep = find(K > 1);

Data.ixSignal_Final = Data.ixSignal_comp(ixKeep);           % the final set of nodes IDs in the signal matrix
Data.ixSignal_Leaves = Data.ixSignal_comp(ixLeaves);        % the leaf nodes
Data.Asignal_final = Data.Asignal_comp(ixKeep,ixKeep);      % the signal network's adjacency matrices

%% save all the things you'll need
save(['Rejected_' fname],'Rejection','Data','pars','optionsModel','optionsReject')


%% look at stuff

% (1) Data eigenvalues against estimated upper and lower bounds
Edata = eig(B);
Retained = Edata > Data.EigEst(1);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 3]); 
line([Data.EigEst(1), Data.EigEst(1)],[0 0.02],'Color',[0.6 0.6 0.8],'Linewidth',1); hold on % upper bound
line([Data.NEigEst(1), Data.NEigEst(1)],[0 0.02],'Color',[0.6 0.6 0.8],'Linewidth',1); hold on % lower bound
line([0,0],[0 0.02],'Color',[0.8 0.8 0],'Linewidth',1); % modularity matrix
plot(Edata,zeros(numel(Edata),1)+0.01,'x','MarkerSize',5,'Color',[0.7 0.7 0.7])
plot(Edata(Retained),zeros(sum(Retained),1)+0.01,'x','MarkerSize',5,'Color',[0.8 0.4 0.4])

lim = max(max(abs(Edata)),Data.EigEst(1)) + 2;
set(gca,'XLim',[-lim lim])
% ylabel('P(max model eigenvalue)')
xlabel('Eigenvalues')

% (2) projection of nodes into first one or two dimensions
nodeIDs = 1:size(Data.A,1);
ixReject = [Rejection.ixNoise; Data.ixSignal_Leaves];

if Data.Dn == 1
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 2]); 
    line([0,0],[0 0.02],'Color',[0 0 0],'Linewidth',1); % modularity matrix
    plot(Data.Dspace(ixReject),zeros(numel(ixReject),1)+0.01,'o','MarkerSize',3,'Color',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]); hold on
    plot(Data.Dspace(Data.ixSignal_Final),zeros(numel(Data.ixSignal_Final),1)+0.01,'o','MarkerSize',3,'Color',[0.8 0.4 0.4],'MarkerFaceColor',[0.8 0.4 0.4])

else
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 6]);
    if ~isempty(ixReject) 
        plot(Data.Dspace(ixReject,1),Data.Dspace(ixReject,2),...
        'o','MarkerSize',3,'Color',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]); hold on
    end
    plot(Data.Dspace(Data.ixSignal_Final,1),Data.Dspace(Data.ixSignal_Final,2),...
        'o','MarkerSize',3,'Color',[0.8 0.4 0.4],'MarkerFaceColor',[0.8 0.4 0.4]); hold on
    xlabel('Dimension 1')
    ylabel('Dimension 2')
end


%% for clustering, see ClusterNetwork.m in /Scripts_for_Paper/





