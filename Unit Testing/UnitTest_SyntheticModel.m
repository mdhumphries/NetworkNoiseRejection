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
alphaList = [-1.5,-1,-0.5,0,0.5,1,1.5];

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
C = cell(numel(NList)*numel(PwithinList)*numel(alphaList),5);
ctr = 0;
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
            ctr = ctr + 1;
            % create weights
            pars = {A,NList{iN},S,alphaList(iA),P};
            [blnExpected,W,blnWithin,blnBetween] = doUnitTest('weight_edges',pars);
            % SANITY CHECKING: check extremes of weights...
            within_wgts = [min(W(blnWithin)) max(W(blnWithin))];
            between_wgts = [min(W(blnBetween)) max(W(blnBetween))];
            
            C(ctr,:) = {NList{iN}, alphaList(iA), within_wgts, between_wgts, sum(sum(W))};
%             if ~blnExpected
%                 disp(['N: ' num2str(NList{iN}) ' Alpha: ' num2str(alphaList(iA)) ...
%                         ' W(within):' num2str(within_wgts) ' W(between):' num2str(between_wgts) ' . Total: ' num2str(sum(sum(W)))])
%             end
        end
        
    end
end

Tresults = cell2table(C,'VariableNames',{'N','Alpha','W_Within','W_Between','Total_W'});