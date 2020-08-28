function [p] = oneTailPermutTestPvalue(sampleObs, loc2Test, tail2Test, nSims)
% One-sample, one-tailed permutation test from Good, "Permutation, Parametric,  
% and Bootstrap Tests of Hypotheses", 2005, p. 34. Testing for the hypothesis 
% that the (unknown) location parameter of the distribution of a given 
% sample is different, and so inconsistent with a given location parameter 
% to test. Analogous to a one-sample, one-tailed t-test for location, but 
% without the need to assume that the data is distributed as a Gaussian.
%
% INPUTS:
% sampleObs: set of observations (those of the unknown location)
% loc2Test: location to test
% tail2Test: if 'right'/'left', the p-value for the location of sample being 
%            larger/smaller than loc2Test is computed. Defaults to 'right'
%            if tail2Test is not recognized.
% nSims: optional number of simulations to run, when the observations in
%        sampleObs is larger than 16, to avoid cluttering RAM. If there are
%        up to 16 observations, an exact p-value is given and nSims is
%        ignored. Otherwise, an approximate p-value is computed based on
%        the nSims simulations.
%
% OUTPUT:
% p: p-value of permutation test
%
% Javier Caballero 13 July 2020



% number of observations
nObs = numel(sampleObs);

% compute deviations
devs = sampleObs - loc2Test;

% shift sample by loc2Test
sampleObs = devs;
loc2Test = 0;

% compute deviations for...
% ... left tail test
if isequal(tail2Test, 'left')
    % sum of negative deviations
    sumDevs = sum(devs(devs < 0));
    
else% ... right tail test by default
    % sum of positive deviations
    sumDevs = sum(devs(devs > 0));
    
end

% if the number of observations in sample is "small"...
if nObs <= 16% ... generate all possible permutations (gets an exact p-value but may get RAM-heavy)
    % all possible binary sequences of [-1 1] of length nObs
    signReassignments = dec2bin(0:2^nObs-1) - '0';
    signReassignments(signReassignments == 0) = -1;% replace 0s by -1s
    
    % sign reassignments applied to the sample
    signReassignments = signReassignments.*repmat(abs(sampleObs), 2^nObs, 1);
    
    % deviations of rearrangements
    devs4Rearrangements = signReassignments - loc2Test;
        
    if isequal(tail2Test, 'left')
        % sum of negative deviations
        sumDevs4Rearrangements = sum(devs4Rearrangements.*(devs4Rearrangements < 0), 2);
        % p-value
        p = sum(sumDevs4Rearrangements <= sumDevs)/(2^nObs);
        
    else
        % sum of positive deviations
        sumDevs4Rearrangements = sum(devs4Rearrangements.*(devs4Rearrangements > 0), 2);
        % p-value
        p = sum(sumDevs4Rearrangements >= sumDevs)/(2^nObs);
        
    end
    
    % report in command window that an exact p-value was computed    
    warning(strcat('The sample to test is SMALL enough (n = ', num2str(nObs),...
        ' <= 16), so the p-value computed is EXACT.'))
    
else% ... use a random reasignment in a given number of simulations
    % pre-allocate
    signReassignments = zeros(nSims, nObs);
    
    % add random sequence of -1 or 1 to matrix of reassignments
    for countSim = 1:nSims
        signReassignments(countSim, :) = randsample([-1 1], nObs, true);
    end
    
    % sign reassignments applied to the sample
    signReassignments = signReassignments.*repmat(abs(sampleObs), nSims, 1);
    
    % deviations of rearrangements
    devs4Rearrangements = signReassignments - loc2Test;
        
    if isequal(tail2Test, 'left')
        % sum of negative deviations
        sumDevs4Rearrangements = sum(devs4Rearrangements.*(devs4Rearrangements < 0), 2);
        % p-value
        p = sum(sumDevs4Rearrangements <= sumDevs)/nSims;
        
    else
        % sum of positive deviations
        sumDevs4Rearrangements = sum(devs4Rearrangements.*(devs4Rearrangements > 0), 2);
        % p-value
        p = sum(sumDevs4Rearrangements >= sumDevs)/nSims;
        
    end
    
    % report in command window that an approximated p-value was computed    
    warning(strcat('The sample to test is LARGE (n = ', num2str(nObs) , ...
        ' > 16), so the p-value computed is Monte Carlo APPROXIMATED.'))
end