% quick script to patch error in saving clustering results...
clearvars

% P(between) = 0.05
fname = 'P_rejection_SyntheticEqual_NoNoise_20180717T123206'; % P(between) = 0.05
run P_Rejection_02_Cluster_NoNoise_synthetic_model

% P(between) = 0.15;
fname = 'P_rejection_SyntheticEqual_NoNoise_20180717T160828'; % P(between) = 0.15
run P_Rejection_02_Cluster_NoNoise_synthetic_model
 
% Unequal groups
fname = 'P_rejection_SyntheticUnequal_NoNoise_20180717T194810';
run P_Rejection_02_Cluster_NoNoise_synthetic_model
