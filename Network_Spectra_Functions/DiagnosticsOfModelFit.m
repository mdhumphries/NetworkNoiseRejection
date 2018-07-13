function diagnostics = DiagnosticsOfModelFit(A,Asamp)

% DIAGNOSTICSOFMODELFIT compares the statistics of network and generated model
% 
% D = DIAGNOSTICSOFMODELFIT(A,M) given n x n weight matrix W of data
% network, and equivalent network M sampled from generative model, returns
% the statistics of various comparisons in struct D, with fields:
%       SAp :  strength distribution of M
%       dS : absolute difference in strength for each node between A and M
%       dSN : difference as fraction of original strength
%       minW :  minimum weight in M
%       maxW : maximum weight in M
%       dmax : absolute difference in maximum weight between A and M
%       MDensity : density of M network (fraction of all entries)
%       dDensity: absolute difference in density between A and M
%
% 13/07/2018: initial version
%
% Mark Humphries
n = size(A,1);
SA = sum(A);
density = sum(A > 0) ./ n^2;

diagnostics.SAp = sum(Asamp);  % degree
diagnostics.minW = min(Asamp); % minimum weight
diagnostics.maxW = max(Asamp);    % maximum weight


% figure; ecdf(kA); hold on; ecdf(kAp); title('Degree distributions of original and permuted network')

diagnostics.dS = abs(SA - diagnostics.SAp);
diagnostics.dSN = 100* diagnostics(iN).dS ./ SA; % difference as fraction of original degree

diagnostics.dmax =  abs(max(A) - diagnostics.maxW);
% figure; ecdf(dKN); title('ECDF of error as proportion of original degree')

diagnostics.MDensity = sum(Asamp > 0) ./ n^2;

diagnostics.dDensity = abs(density - MDensity);

