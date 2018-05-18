function thresh = otsu1D(P)

% OTSU1D Otsu's methof for finding a threshold on bimodality in histograms
% T = OTSU1D(P) given the L-length histogram of discrete probabilities in P,
% iteratively determines the bin number T at which the inter-class variance
% is maximised. T is then the most likely threshold between the two modes
% (if the modes are clearly separated): bins(1:T) are mode 1; bins(T+1:L)
% are mode 2.
% 
% Reference:
% https://en.wikipedia.org/wiki/Otsu%27s_method
%
% Change log:
% 18/5/2018 Initial version
% Mark Humphries 

T = sum(P);
L = numel(P);
if T ~= 1 P = P ./ sum(P); end
T = sum(P);

maximum = 0;
thresh = 1;
w0 = 0; 
for t = 1:L-1    % advance threshold
    w0 = w0 + P(t);     % sum of bins in first mode
    w1 = T - w0;        % sum of bins in second mode
    
    mu0 = sum((1:t) .* P(1:t) / w0);       % current mean of first mode
    mu1 = sum((t+1:L) .* P(t+1:L) / w1);       % current mean of second mode
    
    var_b = w0*w1*(mu0 - mu1).^2;
    
    
    if var_b > maximum
        thresh = t;
        maximum = var_b;
    end
end