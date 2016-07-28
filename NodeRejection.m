function [D,varargout] = NodeRejection(B,Emodel,I,Vmodel,varargin)

% NODEREJECTION separates nodes into "signal" and "noise"
% D = NODEREJECTION(B,E,I,V) splits the nodes in a network into signal 
% and noise components, given: 
%       B: the (nxn) modularity matrix of the network, defined using a null model (e.g Weighted Configuration Model)
%       E: the null-model eigenvalue distribution (n x #repeats of null model) (from e.g. WeightedConfigModel) 
%       I: specified rejection interval (propotion: 0.05, for 95%; 0.01 for
%       99%, and so on); if I is specified as an n-length array {I1,I2,...,In], 
%       then a decompositin will be returned for each I 
%       V: the null model set of eigenvectors (n x n x #repeats of null model; 1 eigenvector per column)
%
% Returns: D, an n-length struct array with fields:
%               .ixSignal: the node indices in the signal component of the
%               network
%               .ixNoise: the node indices in the noise component of the
%               network
%               .N: the number of retained eigenvectors according to I(i) -
%               for use in e.g. estimating the maximum number of communties as N+1
%
% ... = NODEREJECTION(...,Options) passes finds the "noise" nodes by
% weighting X:
%           'none': uses Euclidean distance of projection from the origin
%           'linear': weights projections by the eigenvalues of each
%           eigenvector [Default]
%           'sqrt': weights projections by the square root of the
%           eigenvalues of each eigenvector
%
% Notes: 
% (1) determines the low dimensional space for projecting the network for the specified rejection interval;
% then finds the sampling distribution of null model projections (or norm)
% of each node in that space
% Retains all nodes that exceed the mean+2SEM of their sampling distribution
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
% Mark Humphries 

% sort out options
Options.Weight = 'linear';

N = size(Vmodel,3); 
n = size(Vmodel,1);

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
    [Vpos,ixpos,~] = LowDSpace(B,Emodel,I(i));
    
    %% project data and model
    VmodelW = zeros(n,N);
    switch Options.Weight
        case 'none'
            % Euclidean distance
            Vweighted = Vpos;   % data
            % for each model network, project into the same P dimensions
            % (top P eigenvalues)
            for iN = 1:N   
                VmodelW(:,iN) = sqrt(sum(Vmodel(:,ixpos,iN).^2,2));
            end
        case 'linear'
            egMat = repmat(egs(ixpos)',n,1);
            Vweighted = Vpos .* egMat;  %weight by eigenvalues
            % now do the same for each model repeat: projection, weighted
            % by eigenvalues
            for iN = 1:N
                egMat = repmat(Emodel(ixpos,iN)',n,1);
                VmodelW(:,iN) = sqrt(sum((egMat.*Vmodel(:,ixpos,iN)).^2,2));
            end
        case 'sqrt'
            % weight by square root of eigenvalue: cf Zhang & Newman 2015
            % Phys Rev E
            egMat = repmat((sqrt(egs(ixpos)))',n,1);
            Vweighted = Vpos .* egMat;
            for iN = 1:N
                egMat = repmat((sqrt(Emodel(ixpos,iN)))',n,1);
                VmodelW(:,iN) = sqrt(sum((egMat.*Vmodel(:,ixpos,iN)).^2,2));
            end          
        otherwise
            error('Unknown weighting option')
    end
    
    % do projections
    lengths = sqrt(sum(Vweighted.^2,2));  % length of projection into space
    mModel = mean(VmodelW,2); 
    semModel = std(VmodelW,[],2) / sqrt(N);
    
    % split into signal and noise node sets
    D(i).ixSignal = find(lengths >= mModel + 2.*semModel);  % the retained node
    D(i).ixNoise = find(lengths < mModel + 2.*semModel); % removed nodes
    
%     D(i).ixSignal = find(lengths >= mModel); % + 2.*semModel);  % the retained node
%     D(i).ixNoise = find(lengths < mModel); % + 2.*semModel); % removed nodes

end



