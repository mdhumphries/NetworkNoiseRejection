function S = sample_strength(N,Spar)

% SAMPLE_STRENGTH assign strength sequence to network
% S = SAMPLE_STRENGTH (N,P) creates a stength sequence for a network, given 
% a m-length array N, each entry giving the number of nodes for the part of
% the network with the specified strength distribution.
% S is an m-length struct, with fields:
%       S(i).distribution :         % type of distribution for the ith
%                                       group of nodes
%       S(i).a            :         % First parameter (usually scale) for that distribution
%       S(i).b            :         % Second parameter (usually spread) for that distribution
%
% Returns the n-length array S, giving the strength sequence for the network 
% (where n = sum(N), the total number of nodes)   
%
% Supported distributions (positive support): 
%       'Lognormal': a = mu; b = sd;
%       'Gamma': a = alpha; b = beta;
%       'Poisson': a = lambda (no b);
%       'Binomial': a = N; b = p;
%
% Mark Humphries 7/8/2017

n = sum(N);
m = numel(N); 
S = zeros(n,1);

iN = [0 cumsum(N)];

for iM = 1:m
    % for each case: make column vectors of a 
    a = zeros(N(iM),1) + Spar(iM).a;

    % case call to random with specified distribution
    switch Spar(iM).distribution
        case {'Poisson'}
            % call "random" and assign numbers to array entries.
            S(1+iN(iM):iN(iM+1)) = random(Spar(iM).distribution,a);
            
        case {'Lognormal','Gamma','Binomial'} 
            b = zeros(N(iM),1) + Spar(iM).b;    % two-parameter distributions

            % call "random" and assign numbers to array entries. 
            S(1+iN(iM):iN(iM+1)) = random(Spar(iM).distribution,a,b);
        
        otherwise
            error('Unknown strength distribution specified')
    end
end
