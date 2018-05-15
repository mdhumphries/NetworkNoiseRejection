function [Cs,Cn,ixRetain,CsRetain,CnRetain,D] = deconstructCxy(C,T,varargin) 

% DECONSTRUCTCXY deconstructs a correlation matrix
% [S,N,R,sigR,noiseR,D] = DECONSTRUCTCXY(C,T) deconstructs a matrix Cxy of 
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
%   D: the number of dimensions (eigenvalues above upper bound)
%
%  ... = DECONSTRUCTCXY(...,'weighted') weights the D-dimensional
%  projection by the eigenvalue of each dimension. d = sqrt(sum([eg(i) * V(i)].^2))
%
% NOTES:
%   (1) The correlation/covariance matrix must have zero mean and unit
%   variance. e.g. computed from Z-scored time-series, or explicitly
%   normalised. See MacMahon & Garlaschelli (2015) [their Section II]
%    
%   (2) The Marchenko-Pasteur bounds are taken from Plerou et al 2002 &  MacMahon &
%   Garlaschelli (2015) [their Eq 12]
%
%   (2) Signal + noise decomposition is from MacMahon &
%   Garlaschelli (2015) [their Eqs 13-17]: we obtain the eigenvalues of Cxy, and keep only 
%   those eigenvalues (and corresponding eigenvectors)
%   greater than the upper limit of the Marchenko-Pastur distribution for random
%   matrices. Using these eigenvectors, we then obtain the "signal"
%   correlation matrix C(s) from their outer product. The outer product of
%   all remaining eigenvectors gives us the "noise" correlation matrix C(n)
%
%   (3) Rejecting time-series:  eigenvalues beyond the upper limit indicate 
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
%   R eigenvectors (cf suggestion by Lopes-dos-Santos et al 2011). Option
%   'weighted' weights these projections by the eigenvalue of each
%   dimension (so high variance dimensions account for more of the
%   projection).
%
%   (4) The signal correlation matrix is a further step: we take the
%   sub-matrix of the Correlation matrix for just the retained nodes, 
%   then decompose that into signal and noise components as above.
%
% References:
%   Plerou, V.; Gopikrishnan, P.; Rosenow, B.; Amaral, L. A. N.; Guhr, T. & Stanley, H. E. (2002) 
%   Random matrix approach to cross correlations in financial data. 
%   Phys Rev E, 65, 066126
%
%   Lopes-dos-Santos, V., Conde-Ocazionez, S., Nicolelis, M. A. L., 
%   Ribeiro, S. T. & Tort, A. B. L. (2011) Neuronal assembly detection and cell membership 
%   specification by principal component analysis. PLoS One, 6, e20996
%
%   MacMahon, M. & Garlaschelli, D. (2015) Community Detection for Correlation Matrices 
%   Phys. Rev. X,5, 021006
%
%   14/8/2017: added option to weight projections by the eigenvalue of each
%   dimension
%
% Mark Humphries 7/8/2015

blnWeighted = 0; 
if nargin > 2 && strcmp(varargin{1},'weighted')
    blnWeighted = 1;
end

N = size(C,1);

if N > T
    error('N > T, so correlation matrix is not full rank')
end

% eigenspectra
[V,egs] = eig(C);
egs = diag(egs);
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
if blnWeighted
    egMat = repmat(egs(ixpos)',N,1);
    lengths = sqrt(sum((egMat.*Vpos).^2,2));  % weighted length of projection into space
else
    lengths = sqrt(sum(Vpos.^2,2));  % length of projection into space
end
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

D = numel(ixpos);

% % resulting correlation and noise matrices
% CsRetain = Cs(ixRetain,ixRetain);  % retained signal matrix
% CnRetain = Cn(ixRetain,ixRetain);  % retained noise matrix

