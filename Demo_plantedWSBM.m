%% script to demonstrate planted Weighted Stochastic Block Model
clear all; close all;

% with default parameters
nodes_group=50;
num_groups=2;
mix='low';
frac_periphery=0.5;

options.weight_dist=struct('ingroup', [100,20], 'outgroup', [1,0.02],...
    'periphery_periphery', [1,0.02], 'periphery_core', [10,2]);

    
network = test_noise_rejection_planted_noise(nodes_group,num_groups,mix,frac_periphery,options);