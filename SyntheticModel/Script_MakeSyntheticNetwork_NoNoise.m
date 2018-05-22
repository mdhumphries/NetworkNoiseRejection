%% test script to make the synthetic model without noise
%
% A two-part form of a weighted stochastic block model
% (1) Edges: defined by P-link for class of nodes (planted partition model)
% (2) Weights: generative model using sampled degrees  
% 
% Change log:
%  3/8/2017: initial script to test synthetic model generation 
% 18/5/2018: added estimate of baseline alpha for symmetric distribution of weights, given links
%
% Mark Humphries 

clear all; close all;

addpath ../Network_Analysis_Functions/

% N = [100,100,100];  % size of modules
N = [100,100,100,100];  % size of modules
% N = [100,100];
% N = [200,75,25];  % size of modules

% edge parameters
P.in = 0.1;
P.between = 0.1;

% strength distribution parameters
Spar.distribution =  'Poisson';  % type of distribution
Spar.a = 200;                    % scale: in addition to existing edges
Spar.b = 1;                     % spread

% proportion of degree assigned within module
alpha = 0;    % [-1,1]: -1 = all between; 0 = as per A; 1 = all within

%% check expectations
T = sum(N);  % total nodes
E_in_group = (N-1)/2 .* P.in;      % expected number of unique links in group per node
E_out_group = (T-N)/2 .* P.between; % expected number of unique links outside of group

%% make adjacency matrix
A = wire_edges(N,P);

ixs = 1:sum(N);
Ns = [0 cumsum(N)];

for i = 2:numel(Ns)
    in = ixs > Ns(i-1) & ixs <= Ns(i);
    out = ~in;
    Kwithin(i-1) = mean(sum(A(in,in))) / 2;   % mean degree (unique links)
    Kbetween(i-1) = mean(sum(A(:,in)) - sum(A(in,in))) / 2;
   
end
B = A - expectedA(A);
figure; imagesc(B); colormap(flipud(hot)); colorbar
title('B of A')


%% sample strength distribution according to current rules
S = sample_strength(T,Spar);

%% use Poisson generative model to create Weight matrix

[W,blnWithin,blnBetween] = weight_edges(A,N,S,alpha,P); % calling code from Poisson CM model


% check correct weights
within_wgts = [min(W(blnWithin)) max(W(blnWithin))]
between_wgts = [min(W(blnBetween)) max(W(blnBetween))]


for i = 2:numel(Ns)
    in = ixs > Ns(i-1) & ixs <= Ns(i);
    out = ~in;
    Wbln = W > 0;
    Swithin(i-1) = mean(sum(W(in,in))) / 2;   % mean strength (unique links)
    KWwithin(i-1) =  mean(sum(Wbln(in,in))) / 2; 
    Sbetween(i-1) = mean(sum(W(:,in)) - sum(W(in,in))) / 2;
    KWbetween(i-1) = mean(sum(Wbln(:,in)) - sum(Wbln(in,in))) / 2;
   
end

figure; imagesc(W); colormap(flipud(hot)); colorbar
title('W')

Bw = W - expectedA(W);
figure; imagesc(Bw); colormap(flipud(hot)); colorbar
title('B of W')

