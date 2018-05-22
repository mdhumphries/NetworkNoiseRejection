% unit-test functions for constructing synthetic networks
clear all; close all
addpath('../Network_Analysis_Functions/')
addpath('../SyntheticModel/')

%% fixed parameters
Spar.distribution =  'Poisson';  % type of distribution
Spar.a = 200;                    % scale: in addition to existing edges
Spar.b = 1;                     % spread

%% parameter options for entire model construction
NList = {[100,100],[100 100 100 100],[200,75,25]};
PwithinList = [0.1 0.2 0.5];
P.between = 0.1;
alphaList = [-1,-0.5,0,0.5,1];

%% parameter options for sample strength creation
DistributionList = {'Cauchy','Poisson','Lognormal','Gamma','Binomial'};
a_List = [-5 0 100];
b_List = [-5 0 100];

%% stress test options in sample strength function
for iD = 1:numel(DistributionList)
    for ia = 1:numel(a_List)
        for ib = 1:numel(b_List)
            for iN = 1:numel(NList)
                T = sum(NList{iN});
                thesepars.distribution = DistributionList{iD};
                thesepars.a = a_List(ia);
                thesepars.b = b_List(ia);
                pars = {T,thesepars};
                [blnExpected,S] = doUnitTest('sample_strength',pars);
            end
        end
    end
end

%% construction of synthetic models - each function in turn
for iN = 1:numel(NList)
    % sample strengths
    T = sum(NList{iN});
    pars = {T,Spar};
    [blnExpected,S] = doUnitTest('sample_strength',pars);

    for iP = 1:numel(PwithinList)
        % make adjacency matrix according to block model
        P.in = PwithinList(iP);
        pars = {NList{iN},P};
        [blnExpected,A] = doUnitTest('wire_edges',pars);
        
        for iA = 1:numel(alphaList)
            % create weights
            pars = {A,NList{iN},S,alphaList(iA),P};
            [blnExpected,W] = doUnitTest('weight_edges',pars);
        end
        
        % SANITY CHECKING: check extremes of weights...
    end
end

