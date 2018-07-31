function [D,varargout] = NodeRejection(B,Emodel,I,Vmodel,varargin)

% NODEREJECTION separates nodes into "signal" and "noise"
% D = NODEREJECTION(B,E,I,V) splits the nodes in a network into signal 
% and noise components using the bounds on eigenvalues predicted by a  null model, given: 
%       B: the (nxn) modularity matrix of the network, defined using a null model (e.g Weighted Configuration Model)
%       E: the null-model eigenvalue distribution (n x #repeats of null model) (from e.g. WeightedConfigModel) 
%       I: specified prediction or confidence interval (e.g. I = 0.9 for 90%); if I is specified as an n-length array {I1,I2,...,In], 
%       then a decompositin will be returned for each I 
%       V: the null model set of eigenvectors (n x n x #repeats of null model; 1 eigenvector per column)
%
% Returns: D, an n-length struct array with fields:
%               .ixSignal: the node indices in the signal component of the
%               network
%               .ixNoise: the node indices in the noise component of the
%               network
%               .ixNegative: any nodes indices below the lower bound
%               .N: the number of retained eigenvectors according to I(i) -
%               for use in e.g. estimating the maximum number of communties as N+1
%               .Difference.Raw: a vector of differences between the data
%               values and model-derived threshold values for each node (data-model)
%               .Difference.Norm: as .Raw, but normalised by the model
%               value (-> expressed as a proportion of the model)
%               .NegDiff.Raw: a vector of differences for each node between
%               data and model-derived lower bound (any data < lower bound
%               is interesting)
%               .NegDiff.Norm: as .Raw but normalised as above
%
% ... = NODEREJECTION(...,Options) is a struct that sets analysis options:
%           .Norm: defines the vector norm used to identify noise nodes:
%                'L2': L2-norm AKA the Euclidean distance from the origin in the defined sub-space
%                       [default]
%                'L1': L1-norm AKA the sum of absolute values of the vector
%                'Lmax': L-infinity norm AKA the maximum value
%
%           .Weights: finds the "noise" nodes by weighting X:
%               'linear': weights projections by the eigenvalues of each
%               eigenvector [Default]
%               'none': no weighting by eigenvalues
%               'sqrt': weights projections by the square root of the
%               eigenvalues of each eigenvector
%
%           .Intervals: defines the intervals used to test eigenvalues and
%                   projection lengths: 
%               'PI': prediction intervals on the distribution of values from the null model 
%                       - see PREDICTIONINTERVALSNONP
%               'CI': [Default] confidence intervals on the mean values from the null
%               model (if using this option, set I=0 to just use the mean)
%
%           .Bounds: defines which null model boundary to use:
%               'Upper': creates node projections in space defined by eigenvalues that 
%                        pass the upper-bound of the null model prediction (for e.g. modular structures) [default]
%
%               'Lower': uses the lower-bound of the null model prediction (for e.g. k-partite structures) 
%
% Notes: 
% (1) determines the low dimensional space for projecting the network for the specified rejection interval;
% then finds the sampling distribution of null model projections (or norm)
% of each node in that space
% Retains all nodes that exceed the specified confidence interval of their sampling distribution
%    
% (2) The projections can be weighted according to the eigenvalues (see
% above)
%
% ChangeLog:
% 25/7/2016: initial version
% 26/7/2016: node rejection now based on sampling distribution of the
%               projections
% 28/7/2016: return the dimensionality of the data projection
%               fixed weighted projection bug for the null models
% 29/7/2016: added Norm options; returned difference between model and data  
% 25/10/2016: added checking of lower bound threshold too (assumes there is
% one...)
% 24/7/2017: fixed multiple CI bug; returns estimates of mean and CI for
% all nodes; fixed bug in non-default weighting options
% 25/7/2017: added Prediction Intervals  
% 31/7/2018: added option to project in space defined by lower-bound
%
% Mark Humphries 

addpath('../Helper_Functions/')  % for empty_struct and CIfromSEM

% sort out options
Options.Weight = 'linear';
Options.Norm = 'L2';
Options.Intervals = 'CI';
Options.Bounds = 'Upper';

N = size(Vmodel,3); 
n = size(Vmodel,1);

% update option field values
if nargin > 4
    if isstruct(Options) 
        tempopts = varargin{1}; 
        fnames = fieldnames(tempopts);
        for i = 1:length(fnames)
            Options = setfield(Options,fnames{i},getfield(tempopts,fnames{i}));
        end
    end
end

% compute eigenvalues of modularity matrix
egs = sort(eig(B),'descend');

% get node rejections....
D = emptyStruct({'ixSignal','ixNoise'},[numel(I) 1]);
for i = 1:numel(I)
    % find bounds, and calculate dimensions to retain
    % [Vpos,ixpos,~] = LowDSpace(B,Emodel,I(i));
    
    % choose which projection to use
    switch Options.Bounds
        case 'Upper'
                [Vs,ixVs,~,~] = LowDSpace(B,Emodel,I(i));
                disp('Hello, upper bound here')
        case 'Lower'
                [~,~,~,~,Vs,ixVs] = LowDSpace(B,Emodel,I(i));
                disp('Hello, lower bound here')
        otherwise 
            error('Unknown option for Bounds')
    end
    
    %% project data and model
    nPos = numel(ixVs);
    VmodelW = zeros(n,nPos,N);
    switch Options.Weight
        case 'none'
            % no weighting
            Vweighted = Vs;   % data
            % for each model network, project into the same P dimensions
            % (top P eigenvalues)
 
            VmodelW = Vmodel(:,ixVs,:);

        case 'linear'  % default
            egMat = repmat(egs(ixVs)',n,1);
            Vweighted = Vs .* egMat;  %weight by eigenvalues
            % now do the same for each model repeat: projection, weighted
            % by eigenvalues
            for iN = 1:N
                egMat = repmat(Emodel(ixVs,iN)',n,1);
                VmodelW(:,:,iN) = egMat.*Vmodel(:,ixVs,iN);
            end
        case 'sqrt'
            % weight by square root of eigenvalue: cf Zhang & Newman 2015
            % Phys Rev E
            egMat = repmat((sqrt(egs(ixVs)))',n,1);
            Vweighted = Vs .* egMat;
            for iN = 1:N
                egMat = repmat((sqrt(Emodel(ixVs,iN)))',n,1);
                VmodelW(:,iN) = egMat.*Vmodel(:,ixVs,iN);
            end          
        otherwise
            error('Unknown weighting option')
    end
    
    % norms
    VmodelL = zeros(n,N);
    switch Options.Norm
        case 'L2'  % default
            lengths = sqrt(sum(Vweighted.^2,2));  % L2: Euclidean distance
            for iN = 1:N   
                VmodelL(:,iN) = sqrt(sum(VmodelW(:,:,iN).^2,2));
            end
        case 'L1'
            lengths = sum(abs(Vweighted),2);  % L1: absolute sum
            for iN = 1:N   
                VmodelL(:,iN) = sum(abs(VmodelW(:,:,iN)),2);
            end
        case 'Lmax'
            lengths = max(abs(Vweighted),[],2);  % L infinity: maximum absolute value
            for iN = 1:N   
                VmodelL(:,iN) = max(abs(VmodelW(:,:,iN)),[],2);
            end
       
        otherwise
             error('Unknown Norm')
    end

    
    
    % summarise model projections
    D(i).mModel = mean(VmodelL,2); 

    switch Options.Intervals 
        case 'CI'

            % confidence interval on the mean
            D(i).CIModel = CIfromSEM(std(VmodelL,[],2),ones(size(D(i).mModel,1),1)*N,I(i));

            % differences
            D(i).Difference.Raw = lengths - (D(i).mModel + D(i).CIModel);
            D(i).Difference.Norm = D(i).Difference.Raw ./ (D(i).mModel + D(i).CIModel);
            
            % also store negative projections - only of use if we use CI > 0
            D(i).NegDiff.Raw = lengths - (D(i).mModel - D(i).CIModel);
            D(i).NegDiff.Norm = D(i).NegDiff.Raw ./ (D(i).mModel - D(i).CIModel);

        case 'PI'
            if rem(N,2) == 0
                nIdx = N-1;
            else
                nIdx = N;
            end

            for iN = 1:n
                PIs = PredictionIntervalNonP(VmodelL(iN,1:nIdx));  % take an odd number to (most likely) get integer-spaced prediction intervals
                ix = find(PIs(:,1)/100 == I(i));       % get specified prediction interval 
                if isempty(ix) error('Cannot find specified prediction interval'); end

                D(i).PIModel(iN,:) = PIs(ix,2:3); 
            end
            % differences from upper bound
            D(i).Difference.Raw = lengths - D(i).PIModel(:,2);
            D(i).Difference.Norm = D(i).Difference.Raw ./ D(i).PIModel(:,2);
            
            % also store negative projections below lower bound
            D(i).NegDiff.Raw = lengths - D(i).PIModel(:,1);  % so below lower bound will be negative 
            D(i).NegDiff.Norm = D(i).NegDiff.Raw ./ D(i).PIModel(:,1);

        otherwise
            error('Unknown Interval option')
    end

    % split into signal and noise node sets
    D(i).ixSignal = find(D(i).Difference.Raw > 0);  % the retained node
    D(i).ixNoise = find(D(i).Difference.Raw <= 0); % removed nodes
    D(i).ixNegative = find(D(i).NegDiff.Raw <= 0); % projections below lower bounds
    

end



