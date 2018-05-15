function  [Dspace,ixpos,Dn,mci,varargout] = LowDSpace(B,Emodel,I)

% LOWDSPACE find and return low-dimensional axes for network
% [D,X,N,M] = LOWDSPACE(B,E,I) finds the P low-dimensional axes for the network, given: 
%       B: the (nxn) modularity matrix of the data network, defined using a null model (e.g Weighted Configuration Model)
%       E: the null-model eigenvalue distribution (n x #repeats of model) (from e.g. WeightedConfigModel) 
%       I: specified confidence interval on the maximum eigenvalue
%       (scalar)[e.g. 0.9; enter 0 to just use the mean estimate when using 'CI']
%
%  Returns:
%       D: the n x P matrix of retained eigenvectors (positive eigenvalues)
%       X: the corresponding set of indices into the full eigenvector matrix
%       N: the number of eigenvectors retained
%       M: the model's predicted interval on the maximum eigenvalue. 
%               If 'CI':  [m CI] mean and symmetric confidence interval on the maxium eigenvalue
%               If 'PI':  [L U], lower and upper bound I% on prediction interval
%       
% ... = LOWDSPACE(...,Options) is a struct that sets analysis options:%
%           .Intervals: defines the intervals used to test eigenvalues and
%                   projection lengths: 
%               'PI': prediction intervals on the distribution of values from the null model 
%                       - see PREDICTIONINTERVALSNONP
%               'CI': [Default] confidence intervals on the mean values from the null
%               model (if using this option, set I=0 to just use the mean)
%                       - see CIFROMSEM
%
%
% [...,Dn,Xn,Nn,Mn] = LOWDSPACE(...) are optional outputs for checking for the
% presence of k-partite structure, corresponding to negative data eigenvalues below
% the predicted lower bounds. Each output variable corresponds to its
% counterpart above, but for the negative eigenvalues.
%
%   Notes: a helper function for NodeRejection, but also useful for calling
%   in its own right to just obtain the low-D space
%
%  ChangeLog:   
%   28/7/2016: initial version
%    1/8/2016: input check; implemented maximum eigenvalue rejection
%   25/7/2017: added Prediction Intervals; returned lower bounds too
%   26/7/2017: added options to select interval type
%
% Mark Humphries

Options.Intervals = 'CI';

n = size(B,1);
if size(Emodel,1) ~= n
    error('Eigenvalue matrix should be n (nodes) x N (repeats)')
end

[V,egs] = eig(B);  % eigenspectra of data modularity matrix
egs = diag(egs); % extract vector from diagonal
[egs,ix] = sort(egs,'descend'); % sort eigenvalues into descending order 
V = V(:,ix);  % sort eigenvectors accordingly

% extreme eigenvalues from each model
mx = max(Emodel);       % max (largest positive)
minx = min(Emodel);     % min (largest negative)

switch Options.Intervals 
    case 'CI'
        %% compute mean and CI over largest eigenvalue
        M = mean(mx);
        CIs = CIfromSEM(std(mx),size(Emodel,2),I);
        bnds = M+CIs; % upper confidence interval on maximum eigenvalue for null model

        Mmin = mean(minx);
        CIsMin = CIfromSEM(std(minx),size(Emodel,2),I);
        negbnds = Mmin - CIsMin;  % lower bonds of confidence interval on mean minimum eigenvalue 

        mci = [M CIs];
        varargout{4} = [Mmin CIsMin];
    case 'PI'
        if rem(numel(mx),2) == 0
            n = numel(mx)-1;
        else
            n = numel(mx);
        end
        PIs = PredictionIntervalNonP(mx(1:n));  % take an odd number to (most likely) get integer-spaced prediction intervals
        ix = find(PIs(:,1)/100 == I);       % get specified prediction interval 
        if isempty(ix) error('Cannot find specified prediction interval'); end

        PIsMin = PredictionIntervalNonP(minx(1:n));  % take an odd number to (most likely) get integer-spaced prediction intervals
        ixMin = find(PIsMin(:,1)/100 == I);       % get specified prediction interval 
        if isempty(ix) error('Cannot find specified prediction interval'); end

        % set bounds
        bnds = PIs(ix,3); % upper bounds of positive PI
        negbnds = PIsMin(ixMin,2); % lower bounds of negative PI

        mci = PIs(ix,2:3);
        varargout{4} = PIsMin(ix,2:3);
    otherwise
        error('Unknown Options.Intervals')
end

% option 2: just compute pooled distribution, and return bounds
% prctI = [I/2*100 100-I/2*100]; % rejection interval as symmetric percentile bounds
% bnds = max(prctile(Emodel(:),prctI)); % confidence interval on eigenvalue distribution for null model

%% return dimensions
ixpos = find(egs >= bnds); % eigenvalues exceeding these bounds
ixneg = find(egs <= negbnds);

% return answers
Dn = numel(ixpos);   % number of retained dimensions          
Dspace = V(:,ixpos);  % axes of retained dimensions

varargout{1} = V(:,ixneg);
varargout{2} = ixneg;
varargout{3} = numel(ixneg);

