%% test script to make the synthetic model without noise
%
% A two-part form of a weighted stochastic block model
% (1) Edges: defined by P-link for class of nodes (planted partition model)
% (2) Weights: generative model using sampled degrees  
%
% Mark Humphries 3/8/2017

clear all; close all;

N = [100,100,100];  % size of modules

% edge parameters
P.in = 0.1;
P.between = 0.05;

% degree distribution parameters
Dpar.distribution =         % type of distribution
Dpar.scale = ...            % [in, between]
Dpar.spread = ...           % [in, between]

%% make adjacency matrix
A = wire_edges(N,P);

%% sample degrees according to current rules

D = sample_degrees(N,Dpar);


%% use Poisson generative model to create Weight matrix

W = weight_edges(A,D); % calling code from Poisson CM model
