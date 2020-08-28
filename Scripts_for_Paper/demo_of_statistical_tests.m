%% a demonstration of doing statistical tests on data eigenvalues using the null model's distribution of maximum eigenvalues
% The purpose of this script is to demonstrate ideas from the Supplementary
% Information for future development, on how to extend the spectral estimation approach to include
% statistical testing. 
%
% demonstrates:
% (1) confidence intervals on the expected maximum eigenvalue from the null
% model
% (2) single-sample t-test of each of the data network's largest
% eigenvalues
% (3) single-sample permutation test of each of the data network's largest
% eigenvalues
%
% 28/08/2020 : first version
% Mark Humphries 

clearvars; close all

% path to permutation test function
addpath('Helper_Functions/')

% network to load eigenvalue analysis for
fname = 'LesMis'; 

% parameters
CI = 95;      % which CI? In percentage
nPermutes = 1000;    % how many permutations

%% load data network results, and get maximum eigenvalues
% load 
load(['Results/Rejected_' fname ]);

max_eigs = max(Data.Emodel);

%% confidence interval on the mean value for the maximum eigenvalue

% get CI of mean value
SD = std(max_eigs);
CI_null = CIfromSEM(SD,pars.N,CI/100);  % standard devation of null model, number of data points, size of CI (as a proportion) 

CI_upper_bound = Data.EigEst(1) + CI_null;

% compare to data
B = Data.A - Data.ExpA;  % modularity matrix using chosen null model
Edata = eig(B); Edata = sort(Edata,'descend');

% index of data eigenvalues above upper bound of CI
find(Edata > CI_upper_bound)

%% one-sample t-test

% for each data eigenvalue that exceeds the mean maximum eigenvalue from the null model
% compute a one-sided t-test: we treat the the data eigenvalue as the
% location parameter, and test whether the null model distribution is below
% the location (i.e. a left-sided test)
% t-value = [mean of null model distribution of maximum eigenvalues] - [data eigenvalue] / [SEM of null model distribution of max eigenvalues]

for iE = 1:Data.Dn
    t =  (Data.EigEst(1) - Edata(iE)) / (std(max_eigs) / sqrt(pars.N));
    [~, p_values_t(iE)] = ttest(max_eigs, Edata(iE), 'Tail', 'left');
end

%% permutation test

% for each data eigenvalue that exceeds the E(maximum) from the null model
% compute a one-sided permutation test that the differences betwen the data
% eigenvalue and the null model distribution of maximum eigenvalues exceeds
% zero
for iE = 1:Data.Dn
    p_values_permute(iE) = oneTailPermutTestPvalue(max_eigs, Edata(iE), 'left', nPermutes);
end

% lower limit of P-value (if p-value returns 0, then is less than this
% lower bound):
min_P_value_obtainable = 1/ nPermutes



    