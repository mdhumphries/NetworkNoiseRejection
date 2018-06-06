%% test script to make the synthetic model *with* noise
%
% A two-part form of a weighted stochastic block model
% (1) Edges: defined by P-link for class of nodes (planted partition model)
% (2) Weights: generative model using sampled degrees  
% 
% Change log:
%  3/8/2017: initial script to test synthetic model generation 
% 18/5/2018: added estimate of baseline alpha for symmetric distribution of weights, given links
%  5/6/2018: cloned from NoNoise script to add Noise nodes
% Mark Humphries 

clear all; close all;

addpath ../Network_Analysis_Functions/

% N = [100,100,100];  % size of modules
N = [100,100,100,100];  % size of modules
% N = [100,100];
% N = [200,75,25];  % size of modules

% edge parameters
P.in = 0.2;
P.between = 0.1;

% strength distribution parameters
Spar.distribution =  'Poisson';  % type of distribution
Spar.a = 200;                    % scale: in addition to existing edges
Spar.b = 1;                     % spread

% noise halo parameters
fnoise = 0.5; % proportion of core that is added as noise halo 
P.noise = 0.1;

%% make adjacency matrix
A = wire_edges_noise(N,fnoise,P);


ixs = 1:sum(N);
Ns = [0 cumsum(N)];

for i = 2:numel(Ns)
    in = ixs > Ns(i-1) & ixs <= Ns(i);
    out = ~in;
    Kwithin(i-1) = mean(sum(A(in,in))) / 2;   % mean degree (unique links)
    Kbetween(i-1) = mean(sum(A(:,in)) - sum(A(in,in))) / 2;
    Koutside(i-1) = mean(sum(A(in,out)));
   
end
B = A - expectedA(A);
figure; imagesc(B); colormap(flipud(hot)); colorbar
title('B of A')


%% sample strength distribution according to current rules
T = size(A,1);
S = sample_strength(T,Spar);

%% use Poisson generative model to create Weight matrix
W = weight_edges_noise(A,S); % calling code from Poisson CM model

figure; imagesc(W); colormap(flipud(hot)); colorbar
title('W')

Bw = W - expectedA(W);
figure; imagesc(Bw); colormap(flipud(hot)); colorbar
title('B of W')



