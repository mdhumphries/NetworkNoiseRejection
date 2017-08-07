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
P.between = 0.1;

% stength distribution parameters
Spar.distribution =  'Poisson';  % type of distribution
Spar.a = 10;                    % scale: in addition to existing edges
Spar.b = 1;                     % spread

% proportion of degree assigned within module
alpha = 0.8;    % [0,1): 1 = all within

%% make adjacency matrix
A = wire_edges(N,P);

%% sample strength distribution according to current rules

n = sum(N);
S = sample_strength(n,Spar);


%% use Poisson generative model to create Weight matrix

W = weight_edges(A,N,S,alpha); % calling code from Poisson CM model
