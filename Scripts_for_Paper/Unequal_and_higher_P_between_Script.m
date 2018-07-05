% script to do a sequence of outstanding sets of synthetic tests
% (1) Higher P(between): more dense graphs
%           (a) rejection
%           (b) clustering
% (2) Unequal size groups, using baseline P(between)
%           (a) rejection
%           (b) clustering
%
% Mark Humphries 05/07/2018

%% 1.a. P(between) higher: rejection

Model.N = [100,100,100];  % size of modules
Model.P.between = 0.15;   
Model.Spar.a = 200;                    % scale: in addition to existing edges
Model.P_of_within = [Model.P.between:0.025:0.3];

run P_Rejection_01_NoNoise_synthetic_model

% 1.b cluster results: fname is in workspace
run P_Rejection_02_Cluster_NoNoise_synthetic_model

%% 2.a. Unequal groups
Model.N = [200,75,25];  % size of modules
Model.P.between = 0.05;   
Model.Spar.a = 200;                    % scale: in addition to existing edges
Model.P_of_within = [Model.P.between:0.025:0.2];

run P_Rejection_01_NoNoise_synthetic_model

% 1.b cluster results: fname is in workspace
run P_Rejection_02_Cluster_NoNoise_synthetic_model


