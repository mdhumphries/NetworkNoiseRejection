% script to do a sequence of synthetic model tests
% (1) Starting P(between): sparse graphs
% (2) Higher P(between): more dense graphs
%           (a) rejection
%           (b) clustering
% (3) Unequal size groups, using baseline P(between)
%           (a) rejection
%           (b) clustering
%
% (4) Noise nodes on original sparse graphs
%
% 05/07/2018 : original version for just higher and unequal
% 16/07/2018 : all tests  
% 
% Mark Humphries 
clearvars

nModels = 2;

%% 1. P(between) baseline model
Model.P.between = 0.05;   
Model.P_of_within = [Model.P.between:0.025:0.2];
Model.P_of_within = [Model.P.between,0.2];

Model.N = [100,100,100,100];  % size of modules
Model.Spar.a = 200;                    % scale: in addition to existing edges

% 1.a reject
run P_Rejection_01_NoNoise_synthetic_model

% 1.b cluster results: fname is in workspace
run P_Rejection_02_Cluster_NoNoise_synthetic_model


%% 2. P(between) higher
Model.P.between = 0.15;   
Model.P_of_within = [Model.P.between:0.025:0.3];
% 2.a reject
run P_Rejection_01_NoNoise_synthetic_model

% 2.b cluster results: fname is in workspace
run P_Rejection_02_Cluster_NoNoise_synthetic_model

%% 3. Unequal groups
Model.N = [200,100,75,25];  % size of modules
Model.P.between = 0.05;   
Model.P_of_within = [Model.P.between:0.025:0.2];

% rejection
run P_Rejection_01_NoNoise_synthetic_model

% cluster results: fname is in workspace
run P_Rejection_02_Cluster_NoNoise_synthetic_model

%% 4. adding noise nodes to equal sized groups
nModels = 50;

Model.P.in = 0.2;   % IMPORTANT: chose this carefully...
Model.P.between = 0.05;
Model.N = [100,100,100,100];  % size of modules
Model.F_noise = [0.25 0.5 1];
% note: Model.P_of_noise is determined in the script, given all other
% properties....
run P_Rejection_01_Noise_synthetic_model.m

fname = 'P_rejection_SyntheticEqual_Noise_20180716T123952';
% cluster results: fname is in workspace
run P_Rejection_02_Cluster_Noise_synthetic_model.m


