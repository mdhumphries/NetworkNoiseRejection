function [Cs,Cn,ixRetain,CsRetain,CnRetain] = deconstructCxy(C,T) 

% DECONSTRUCTCXY deconstructs a correlation matrix
% [S,N,R,sigR,noiseR] = DECONSTRUCTCXY(C,T) deconstructs a matrix Cxy of 
% nxn pairwise correlation coefficients, computed from n time-series recorded for T 
% time-steps. The deconstruction assumes that the observed correlation 
% matrix is constructed from a true set of correlations (signal) with
% IID added noise from a symmetric distribution. The function estimates this 
% decomposition two ways: on the full correlation matrix; and on a reduced 
% correlation matrix from which time-series that are 
% unrelated to all the others have been discarded. It returns:
%   S: the matrix of estimated signal correlations 
%   N: the matrix of estimated noise correlations (Cxy = S + N)
%   R: the index of all n time-series that are retained   
%   sigR: the matrix of signal correlations between retained time-series 
%   noiseR : the matrix of noise correlations between retained time-series 
%
% NOTES:
%   (1) Signal + noise decomposition is from MacMahon &
%   Garlaschelli (2013): we obtain the eigenvalues of Cxy, and keep only 
%   those eigenvalues (and corresponding eigenvectors)
%   greater than the upper limit of the Marchenko-Pastur distribution for random
%   matrices. Using these eigenvectors, we then obtain the "signal"
%   correlation matrix C(s) from their outer product. The outer product of
%   all remaining eigenvectors gives us the "noise" correlation matrix C(n)
%
%   (2) Rejecting time-series:  eigenvalues beyond the upper limit indicate 
%   time-series with orthogonal contributions to the correlation matrix beyond those predicted by 
%   the random noise model. The number of these eigenvalues thus indicates the 
%   number of possible groups of correlated time-series. Eigenvalues below the lower limit indicate
%   contributions aligned to contributions already accounted for. The total number E 
%   of eigenvalues remaining within Marchenko-Pastur distribution's upper and
%   lower limits are thus not contributing to any correlation structure
%   beyind the random noise model, and can be discarded.
% 
%   We find the retained R = N -E time-series by retaining the R time-series 
%   with the longest projecting vectors in the vector space defined by the 
%   R eigenvectors (cf suggestion by Lopes-dos-Santos et al 2011). 
%
% References:
%   Lopes-dos-Santos, V., Conde-Ocazionez, S., Nicolelis, M. A. L., 
%   Ribeiro, S. T. & Tort, A. B. L. (2011) Neuronal assembly detection and cell membership 
%   specification by principal component analysis. PLoS One, 6, e20996
%
%   MacMahon, M. & Garlaschelli, D. (2013) Community detection for correlation
%   matrices. arxiv 1311.1924
%
%
% Mark Humphries 7/8/2015

N = size(C,1);

if N > T
    error('N > T, so correlation matrix is not full rank')
end

% eigenspectra
[V,egs] = eig(C,'vector');
[egs,~] = sort(egs,'descend'); 

% Marchenko-Pastur distribution bounds
MP_lambda_pos = (1+sqrt(N/T)).^2;
MP_lambda_neg = (1-sqrt(N/T)).^2;

ixpos = find(egs > MP_lambda_pos);
ixneg = find(egs < MP_lambda_neg);  % use for rejection
Tensemble = numel(ixpos) + numel(ixneg);  % total number of time-series

% how many to test?    
Vpos = V(:,ixpos);      % corresponding set of eigenvectors

% reconstruct the "signal" correlation matrix
Cs = zeros(N);
for i=1:numel(ixpos)
    Cs = Cs + egs(ixpos(i)) * V(:,ixpos(i))*V(:,ixpos(i))';   
end
Cn = C - Cs; % noise matrix

% rejection of outliers
lengths = sqrt(sum(Vpos.^2,2));  % length of projection into space
[srt,I] = sort(lengths,'descend');
ixRetain = sort(I(1:Tensemble),'ascend');  % the T retained nodes, in ID order
ixRemove = sort(I(Tensemble+1:end),'ascend'); % removed nodes, in ID order

%%% Do deconstruction on retained Cxy
CRetain = C(ixRetain,ixRetain);
[V,D] = eig(CRetain);
egs =diag(D);
ixpos = find(egs > MP_lambda_pos);
CsRetain = zeros(numel(ixRetain));
for i=1:numel(ixpos)
    CsRetain = CsRetain + egs(ixpos(i)) * V(:,ixpos(i))*V(:,ixpos(i))';   
end
CnRetain = CRetain - CsRetain; % noise matrix


% % resulting correlation and noise matrices
% CsRetain = Cs(ixRetain,ixRetain);  % retained signal matrix
% CnRetain = Cn(ixRetain,ixRetain);  % retained noise matrix

