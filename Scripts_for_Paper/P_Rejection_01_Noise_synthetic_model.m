%% script to assess performance of group detection with noise halo
%
% Here in scrip 01: create networks and detect with spectral rejection
%
% Change log:
% 6/6/2018: initial version
%
% Mark Humphries

clear all; close all;
addpath('../SyntheticModel/');
addpath('../Network_Spectra_Functions/');
fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';

%% fixed parameters of synthetic model
% Model.N = [200,75,25];  % size of modules
Model.N = [100,100,100];  % size of modules


Model.P.in = 0.2;   % IMPORTANT: chose this carefully...
Model.P.between = 0.05;  % has to be same as for No-noise network...   

% stength distribution parameters
Model.Spar.distribution =  'Poisson';  % type of distribution
Model.Spar.a = 200;                    % scale: in addition to existing edges
Model.Spar.b = 1;                     % spread

nBatch = 50;

%% range of parameters
Model.F_noise = [0.25 0.5 1];
% 3 qualitative cases of P(noise)
Model.P_of_noise = [Model.P.between/2, Model.P.between, Model.P.between + (Model.P.in-Model.P.between)/2, Model.P.in, Model.P.in + (Model.P.in-Model.P.between)/2];

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
P.in = Model.P.in;

for iP = 1:numel(Model.P_of_noise)
    % assign parameters   
    P.noise = Model.P_of_noise(iP);
    for iF = 1:numel(Model.F_noise)
        fnoise = Model.F_noise(iF);
        
        for iB = 1:nBatch
            disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_noise)) '; F:' num2str(iF) '/' num2str(numel(Model.F_noise)) '; batch: ' num2str(iB)])
            tic
            %% make model
            A = wire_edges_noise(Model.N,fnoise,P);  % make adjacency matrix
            T = size(A,1);
            S = sample_strength(T,Model.Spar); % sample strength distribution according to current rules            
            W = weight_edges_noise(A,S); % use Poisson generative model to create Weight matrix

            %% do spectral rejection on model
            [Data,Rejection,Control] = reject_the_noise(W,rejectionpars,optionsModel,optionsReject);        

            % store relevant network information for clustering
            Network(iP,iF,iB).W = W;
            Network(iP,iF,iB).ExpW = Data.ExpA;
            Network(iP,iF,iB).WsignalFinal = Data.Asignal_final;
            Network(iP,iF,iB).ixFinal = Data.ixSignal_Final;
            
            % RESULTS:
            % number of groups recovered by spectra vs same count from other approaches

            Results.SpectraWCM.Groups(iP,iF,iB) = Data.Dn+1;
            Results.SpectraConfig.Groups(iP,iF,iB) = Control.Dn+1;
            Results.PosEigWCM.Groups(iP,iF,iB) = Data.PosDn+1;
            Results.PosEigConfig.Groups(iP,iF,iB) = Control.PosDn+1;
            
            % proportion of nodes assigned correctly...
            core = 1:n;           % indices of core nodes
            halo = n+1:T;         % indices of halo nodes
            
            % (1) True positives
            Results.SpectraWCM.TruePos(iP,iF,iB) = numel(intersect(Rejection.ixSignal,core));
            % (2) False positives
            Results.SpectraWCM.FalsePos(iP,iF,iB) = numel(intersect(Rejection.ixSignal,halo));
            % (3) True negatives
            Results.SpectraWCM.TrueNeg(iP,iF,iB) = numel(intersect(Rejection.ixNoise,halo));
            % (4) False negatives
            Results.SpectraWCM.FalseNeg(iP,iF,iB) = numel(intersect(Rejection.ixNoise,core));
            
            % Rates:
            Results.SpectraWCM.Sensitivity(iP,iF,iB) = Results.SpectraWCM.TruePos(iP,iF,iB) ./( Results.SpectraWCM.TruePos(iP,iF,iB) + Results.SpectraWCM.FalseNeg(iP,iF,iB)); % proportion of accepted  
            Results.SpectraWCM.Specificity(iP,iF,iB) = Results.SpectraWCM.TrueNeg(iP,iF,iB) ./( Results.SpectraWCM.TrueNeg(iP,iF,iB) + Results.SpectraWCM.FalsePos(iP,iF,iB));  % proportion of rejected
            
            Results.Time(iP,iB) = toc;
        end
    end
end

%% process results
Results.ProportionModular.SpectraWCM = squeeze(sum(Results.SpectraWCM.Groups > 1,3));
Results.ProportionModular.SpectraConfig = squeeze(sum(Results.SpectraConfig.Groups > 1,3));
Results.ProportionModular.PosEigWCM = squeeze(sum(Results.PosEigWCM.Groups > 1,3));
Results.ProportionModular.PosEigConfig = squeeze(sum(Results.PosEigConfig.Groups > 1,3));

    
%% save
fname = datestr(now,30);

if all(Model.N == Model.N(1))
    strName = 'Equal';
else
    strName = 'Unequal';
end

save([fpath 'P_rejection_Synthetic' strName '_Noise_' fname],'Results','Network','Model','rejectionpars','optionsModel','optionsReject','-v7.3')


