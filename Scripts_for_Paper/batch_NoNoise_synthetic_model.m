%% script to assess performance on equal-sized groups 
% 
% Change log:
% 7/8/2017: initial script with synthetic model and clustering
% 4/6/2018: split into spectral analysis, and clustering
%
% Mark Humphries

clear all; close all;
addpath('../SyntheticModel/');
addpath('../Network_Spectra_Functions/');

%% fixed parameters of synthetic model
% Model.N = [200,75,25];  % size of modules
Model.N = [100,100,100];  % size of modules

Model.P.between = 0.05;   

% stength distribution parameters
Model.Spar.distribution =  'Poisson';  % type of distribution
Model.Spar.a = 50;                    % scale: in addition to existing edges
Model.Spar.b = 1;                     % spread

nBatch = 2;

%% range of parameters
Model.P_of_within = [0.05:0.025:0.2];
% Model.alpha_range = [0.5:0.025:0.6 0.65:0.05:0.8];
Model.alpha_range = [-1 -0.5 -0.1 0 0.1 0.5 1];

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

% later: loop here to resample degree sequence
Model.S = sample_strength(n,Model.Spar); % sample strength distribution according to current rules

for iP = 1:numel(Model.P_of_within)
    for iA = 1:numel(Model.alpha_range)
        % assign parameters   
        P.in = Model.P_of_within(iP);
        alpha = Model.alpha_range(iA);
        
        for iB = 1:nBatch
            disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_within)) '; alpha ' num2str(iA) '/' num2str(numel(Model.alpha_range)) '; batch: ' num2str(iB)])
            tic
            %% make model
            A = wire_edges(Model.N,P);  % make adjacency matrix
            W = weight_edges(A,Model.N,Model.S,alpha,P); % use Poisson generative model to create Weight matrix

            %% do spectral rejection on model
            [Data,Rejection,Control] = reject_the_noise(W,rejectionpars,optionsModel,optionsReject);        
            
            % store relevant network information for clustering
            Network(iP,iA,iB).A = Data.A;
            Network(iP,iA,iB).ExpA = Data.ExpA;
            Network(iP,iA,iB).Asignal_final = Data.Asignal_final;
            Network(iP,iA,iB).ixSignal_Final = Data.ixSignal_Final;
            % RESULTS:
            % number of groups recovered by spectra vs same count from other
            % approaches
            
            Results.SpectraWCMGroups(iP,iA,iB) = Data.Dn+1;
            Results.SpectraConfigGroups(iP,iA,iB) = Control.Dn+1;
            Results.PosEigWCMGroups(iP,iA,iB) = Data.PosDn+1;
            Results.PosEigConfigGroups(iP,iA,iB) = Control.PosDn+1;
            
            Results.Time(iP,iA,iB) = toc;
        end
    end
end

%% save
fname = datestr(now,30);

if all(Model.N == Model.N(1))
    strName = 'Equal';
else
    strName = 'Unequal';
end

save(['../Results/Synthetic' strName '_NoNoise_' fname],'Results','Network','Model','group_membership','rejectionpars','optionsModel','optionsReject')


