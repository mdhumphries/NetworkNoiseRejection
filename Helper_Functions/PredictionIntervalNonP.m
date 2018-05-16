function PIs = PredictionIntervalNonP(X)

% PREDICTIONINTERVALNONP non-parametric prediction interval
% PI = PREDICTIONINTERVALNONP(X) given an n-length array of sampled data X, 
% computes the non-parameteric prediction intervals PI.
%
% PI: n/2 x 3 array, one row per possible PI. Each row is:
%   [PI% L U] : the PI% level, and the lower (L) and upper (U) bounds 
%
% Reference: 
% Preston, S. Teaching Prediction Intervals Journal of Statistics Education, 2000, 8, 3
%
% Mark Humphries 25/7/2017

n = numel(X);
cutoff = floor(n/2);

Xorder = sort(X,'ascend');

PIs = zeros(cutoff,3);
for j = 1:cutoff
    PI = (1-2*j/(n+1)) * 100;
    PIs(j,:) = [PI Xorder(j) Xorder(n+1-j)];
end
