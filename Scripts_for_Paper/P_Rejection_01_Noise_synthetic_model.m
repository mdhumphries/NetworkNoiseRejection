%% script to assess performance of group detection with noise halo
%
% Here in scrip 01: create networks and detect with spectral rejection
%
% Change log:
% 6/6/2018: initial version
%
% Mark Humphries

% clear all; close all;
clearvars Network Results

addpath('../SyntheticModel/');
addpath('../Network_Spectra_Functions/');
fpath = 'C:/Users/lpzmdh/Dropbox/Analyses/Networks/SyntheticModel_Rejection_Results/';

%% parameters to potentially explore - comment out to run batch script

% Model.N = [200,75,25];  % size of modules
% Model.N = [100,100,100];  % size of modules
% Model.Spar.a = 200;                    % scale: in addition to existing edges
% 
% Model.P.in = 0.2;   % IMPORTANT: chose this carefully...
% Model.P.between = 0.05;  % has to be same as for No-noise network...   

%% range of parameters
% Model.F_noise = [0.25 0.5 1];
% 5 qualitative cases of P(noise)
Model.P_of_noise = [Model.P.between/2, Model.P.between, Model.P.between + (Model.P.in-Model.P.between)/2, Model.P.in, Model.P.in + (Model.P.in-Model.P.between)/2];

%% fixed parameters of synthetic model
% stength distribution parameters
Model.Spar.distribution =  'Poisson';  % type of distribution
Model.Spar.b = 1;                     % spread

% nModels = 50;

%% rejection and clustering parameters
rejectionpars.N = 100;           % repeats of permutation
rejectionpars.I = 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
rejectionpars.Model = 'Link';   % or 'WCM' . % which null model
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
        
        for iB = 1:nModels
            disp(['P: ' num2str(iP) '/' num2str(numel(Model.P_of_noise)) '; F:' num2str(iF) '/' num2str(numel(Model.F_noise)) '; batch: ' num2str(iB)])
            tic
            %% make model
            A = wire_edges_noise(Model.N,fnoise,P);  % make adjacency matrix
            T = size(A,1);
            S = sample_strength(T,Model.Spar); % sample strength distribution according to current rules            
            W = weight_edges_noise(A,S); % use Poisson generative model to create Weight matrix

            %% do spectral rejection on model
            [Data,Rejection,Control,ControlRejection] = reject_the_noise(W,rejectionpars,optionsModel,optionsReject);        

            % store relevant network information for clustering
            Network(iP,iF,iB).W = W;
            Network(iP,iF,iB).ExpW = Data.ExpA;
            Network(iP,iF,iB).ixFinal_Sparse = Data.ixSignal_Final;        % index into W to recover final signal matrix from Sparse WCM
            Network(iP,iF,iB).ixFinal_Full = Control.ixSignal_Final;        % index into W to recover final signal matrix from Full WCM 
            
            % RESULTS:
            % number of groups recovered by spectra vs same count from other approaches

            Results.SpectraSparseWCM.Groups(iP,iF,iB) = Data.Dn+1;
            Results.SpectraFullWCM.Groups(iP,iF,iB) = Control.Dn+1;
            Results.PosEigSparseWCM.Groups(iP,iF,iB) = Data.PosDn+1;
            Results.PosEigFullWCM.Groups(iP,iF,iB) = Control.PosDn+1;
            
            % proportion of nodes assigned correctly...
            core = 1:n;           % indices of core nodes
            halo = n+1:T;         % indices of halo nodes
            
            % (1) True positives
            Results.SpectraSparseWCM.TruePos(iP,iF,iB) = numel(intersect(Rejection.ixSignal,core));
            Results.SpectraFullWCM.TruePos(iP,iF,iB) = numel(intersect(ControlRejection.ixSignal,core));
            
            % (2) False positives
            Results.SpectraSparseWCM.FalsePos(iP,iF,iB) = numel(intersect(Rejection.ixSignal,halo));
            Results.SpectraFullWCM.FalsePos(iP,iF,iB) = numel(intersect(ControlRejection.ixSignal,halo));
            
            % (3) True negatives
            Results.SpectraSparseWCM.TrueNeg(iP,iF,iB) = numel(intersect(Rejection.ixNoise,halo));
            Results.SpectraFullWCM.TrueNeg(iP,iF,iB) = numel(intersect(ControlRejection.ixNoise,halo));
            
            % (4) False negatives
            Results.SpectraSparseWCM.FalseNeg(iP,iF,iB) = numel(intersect(Rejection.ixNoise,core));
            Results.SpectraFullWCM.FalseNeg(iP,iF,iB) = numel(intersect(ControlRejection.ixNoise,core));
            
            % Rates:
            Results.SpectraSparseWCM.Sensitivity(iP,iF,iB) = Results.SpectraSparseWCM.TruePos(iP,iF,iB) ./( Results.SpectraSparseWCM.TruePos(iP,iF,iB) + Results.SpectraSparseWCM.FalseNeg(iP,iF,iB)); % proportion of accepted  
            Results.SpectraSparseWCM.Specificity(iP,iF,iB) = Results.SpectraSparseWCM.TrueNeg(iP,iF,iB) ./( Results.SpectraSparseWCM.TrueNeg(iP,iF,iB) + Results.SpectraSparseWCM.FalsePos(iP,iF,iB));  % proportion of rejected
            Results.SpectraFullWCM.Sensitivity(iP,iF,iB) = Results.SpectraFullWCM.TruePos(iP,iF,iB) ./( Results.SpectraFullWCM.TruePos(iP,iF,iB) + Results.SpectraFullWCM.FalseNeg(iP,iF,iB)); % proportion of accepted  
            Results.SpectraFullWCM.Specificity(iP,iF,iB) = Results.SpectraFullWCM.TrueNeg(iP,iF,iB) ./( Results.SpectraFullWCM.TrueNeg(iP,iF,iB) + Results.SpectraFullWCM.FalsePos(iP,iF,iB));  % proportion of rejected

            
            Results.Time(iP,iB) = toc;
        end
    end
end

%% process results
Results.ProportionModular.SpectraSparseWCM = squeeze(sum(Results.SpectraSparseWCM.Groups > 1,3));
Results.ProportionModular.SpectraFullWCM = squeeze(sum(Results.SpectraFullWCM.Groups > 1,3));
Results.ProportionModular.PosEigSparseWCM = squeeze(sum(Results.PosEigSparseWCM.Groups > 1,3));
Results.ProportionModular.PosEigFullWCM = squeeze(sum(Results.PosEigFullWCM.Groups > 1,3));
    
%% save
strDate = datestr(now,30);

if all(Model.N == Model.N(1))
    strName = 'Equal';
else
    strName = 'Unequal';
end

fname = ['P_rejection_Synthetic' strName '_Noise_' strDate];

save([fpath fname],'Results','Network','Model','rejectionpars','optionsModel','optionsReject','-v7.3')


