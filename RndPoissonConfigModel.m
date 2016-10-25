function [E,varargout] = RndPoissonConfigModel(A,N,varargin)

% E = RNDPOISSCONFIGMODEL(A,N) takes the weighted, undirected (n x n) adjacency matrix A and creates N
% random realisations of the modularity matrix B by randomly generating a 
% null model approximating the weighted configuration model P using a Poisson process.
% Returns E, an nxN matrix of all n eigenvalues for all N random modularity
% matrices.
%
% NOTE: if the network has real-valued weights, then the weights are
% converted into multi-edges (i.e. integer weights) by first scaling with
% some factor C; default is 100. This is most appropriate when the weights are e.g.
% correlation values in [0,1]. To omit, or customise, set optional
% conversion factor as outlined below.
%
% ... = RNDPOISSCONFIGMODEL(..,C,OPTIONS) sets optional settings:
%       C: sets the conversion factor C; i.e. the amount by which the 
%           weighted adjacency matrix is scaled to get integer weights.
%           C = 'all' sets the conversion factor large enough that the minimum weight
%           is converted to 1. Set C = [] to omit.
%
%       OPTIONS: struct of options:
%           .Expected = {0,1}: if specified (=1), uses the ensemble of generated
%           configuration models to estimate the expected model. Useful for
%           non-standard configuration models. [default = 0]
%           .NoLoops = {0,1}: if specified (=1), prevents self-loops in the
%           generated random models [default = 1]
%   
% [..,D,V,X] = RNDPOISSCONFIGMODEL(...) returns:
%           D: a struct, containing diagnostic measurements of the accuracy of the null model 
%           for each of the N repeats, with fields
%           D(i).sAp = strength distribution of the ith repeat
%           D(i).dS = absolute difference between data and ith model strength distributions 
%           D(i).dSN = absolute difference, normalised per node to its strength
%               in the data (i.e. to measure the error relative to magnitude)
%           V: an nxnxN matrix, containing all of the nxn eigenvector matrices of the N repeats            
%           X: the nxn matrix for the expected configuration model: only
%           returned if Options.Expected = 1;
%
% Notes: 
% (1) assumes A is connected;
% (2) To Do: Add negative weights as separate option: can split into (+) and (-) groups, and assign 
% links, then weights	
%
% ChangeLog:
% 31/08/2016: changed from weighted configuration model to random
% poissonian generator model. Modified from WEIGHTEDCONFIGMODEL, Mark 
% Humphries's script 26/07/2016 [SM]
% 25/10/16: added options from WEIGHTEDCONFIGMODEL  [MH]
% 25/10/16: updated Poisson model to be equivalent to multinomial draw from
% WCM model [MH]
%
% Silvia Maggi & Mark Humphries

n = size(A,1);

sA = sum(A);  % original strength distribution
kA = sum(A>0);  % original degree distribution

% catch errors
sAc = sum(A'); 
if sA ~= sAc
	warning('Directed networks will be converted to undirected networks');
	A = (A + A') / 2; % convert to undirected
	sA = sum(A); 
	kA = sum(A>0);
end

if any(A < 0)
	error('WCM only defined for positive weights - for now')
end

% update option field values
Options.Expected = 0;
Options.NoLoops = 1;

if nargin >= 4
    if isstruct(Options) 
        tempopts = varargin{2}; 
        fnames = fieldnames(tempopts);
        for i = 1:length(fnames)
            Options = setfield(Options,fnames{i},getfield(tempopts,fnames{i}));
        end
    end
end
% quantisation steps
if nargin >= 3
    conversion = varargin{1};
    if isempty(conversion) conversion = 1; end
    if strfind(conversion,'all')
        % the scale so that minimum non-zero weight is 1
        conversion = 1./minW;
    end
else
    conversion = 100; % into integer number of edges
end

% check if weights are already integers
if ~any(rem(A(:),1))  % then is integers for all weights
    conversion = 1;
end

% convert network into multi-edge version
A_int = round(A*conversion); 

% weighted configuration model [expectation]
P = expectedA(A);  

% initialise structs to full storage size, allowing parfor to slice
% appropriately
Pstar = emptyStruct({'Egs','A','V'},[N,1]);

fields = {'conversion','sAp','minW','maxW','dS','dSN','dmax'};
diagnostics = emptyStruct(fields, [N,1]);

% detect parallel toolbox, and enable if present
blnParallel = license('test','Distrib_Computing_Toolbox');

if blnParallel
    nCores = feature('numCores');
    if isempty(gcp('nocreate'))
        parpool('local',nCores);  % run on all
    end
end


parfor iN = 1:N
% for iN = 1:N
    %% Step 1: create links
    K = sum(kA);  % total number of links
   
    pnode = expectedA(A>0); % probability of link between each pair of nodes
    if Options.NoLoops
        Aperm = real(rand(n,n) < triu(pnode,1));  % don't include diagonal: no self-loops (slightly underestimates degree)
    else
        Aperm = real(rand(n,n) < triu(pnode,0));  % include diagonal: allow self-loops
    end    
    
    %% Step 2: make weights
    S = sum(sA);  % total weights
    sAint = sum(A_int); % integer strength
    Sint = sum(sAint); % integer total strength

    if S ~= K  % then is weighted network  
        % disp('Weighted network')
        % linear indices of linked pairs for upper triangular matrix
        ixpairs = find(triu(Aperm,0)>0); 
        
        % get as (i,j)
        [irow,jcol] = ind2sub([n,n],ixpairs);
        % linear indices of linked pairs for lower triangular matrix
        % (symmetric respect to the diagonal)
        linearInd = sub2ind([n,n], jcol, irow);
        % get P(link): 
        Plink = sA(irow) .* sA(jcol);
        Plink = Plink ./ sum(Plink); % P(link is placed between each pair)
        
        nLinks = round(Sint/2) - numel(ixpairs);  % total links - [already placed]

        lambda = nLinks .* Plink; % expected number of links
        
        Nlink = poissrnd(lambda);  % Poisson random number of links made
        
        % keyboard
        
        Aperm(ixpairs) = Aperm(ixpairs) + Nlink'; % add to existing links
        Aperm = Aperm ./ conversion;  % convert back
        Aperm = Aperm + Aperm';  % make symmetric
        
%         
%         
%         Plink = poissrnd(sA(irow) .* sA(jcol));
%         rndmat(ixpairs) = Plink;
%         rndmat(linearInd) = Plink;
%         Atemp = sum(sum(A))/sum(sum(rndmat))*rndmat; % normalisation
%         % kkk = rndmat.*(ones(size(A,1))-eye(size(A,1))); % remove diagonal
%         % Aperm = sum(sum(A))/sum(sum(kkk))*kkk; % normalisation
        
    end
    
    %% diagnostics: how far does random model depart?
    diagnostics(iN).sAp = sum(Aperm);  % degree
    diagnostics(iN).minW = min(Aperm); % minimum weight
    diagnostics(iN).maxW = max(Aperm);    % maximum weight
    
    % figure; ecdf(sA); hold on; ecdf(sAp); title('Degree distributions of original and permuted network')
    
    diagnostics(iN).dS = abs(sA - diagnostics(iN).sAp);
    diagnostics(iN).dSN = 100* diagnostics(iN).dS ./ sA; % difference as fraction of original degree
 
    diagnostics(iN).dmax =  max(A) - diagnostics(iN).maxW;
    % figure; ecdf(dSN); title('ECDF of error as proportion of original degree')
    
    if Options.Expected 
        Pstar(iN).A = Aperm;  % store permuted network
    else
        %% go ahead and get eigenvalues
        % P is null model for A, assuming A = P + noise
        % B* = P* - P
        [Pstar(iN).V,Pstar(iN).Egs] = eig(Aperm - P,'vector');
        [Pstar(iN).Egs,ix] = sort(Pstar(iN).Egs,'descend'); % ensure eigenvalues are sorted in order
        Pstar(iN).V = Pstar(iN).V(:,ix); % also sort eigenvectors

    end
end

if Options.Expected 
    % generate expected matrix
    Aall = zeros(n);
    for iN = 1:N
      Aall = Aall + Pstar(iN).A;
    end
    ExpWCM = Aall ./ N;
    varargout{3} = ExpWCM;
    % now compute eigenvalues
    for iN = 1:N
        [Pstar(iN).V,Pstar(iN).Egs] = eig(Pstar(iN).A - ExpWCM,'vector');
        [Pstar(iN).Egs,ix] = sort(Pstar(iN).Egs,'descend'); % ensure eigenvalues are sorted in order
        Pstar(iN).V = Pstar(iN).V(:,ix); % also sort eigenvectors
    end
end

% now collapse all eigenvalues and vectors into matrix
V = zeros(n,n,N);
E = zeros(n,N);
A = zeros(n,n,N);
for iN = 1:N
    E(:,iN) = Pstar(iN).Egs;
    V(:,:,iN) = Pstar(iN).V;
    A(:,:,iN) = Pstar(iN).A;
end

varargout{1} = diagnostics;
varargout{2} = V;
% 3 is assigned to the expected WCM above
varargout{4} = A;