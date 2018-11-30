function [E,varargout] = poissonSparseWCM(A,N,varargin)

% E = POISSONSPARSEWCM(A,N) takes the weighted, undirected (n x n) adjacency matrix A and creates N
% random realisations of the modularity matrix B by randomly generating a 
% null model approximating the *sparse* weighted configuration model P using a Poisson process.
%
% The sparse WCM has a two-step approach:
% (1) Assign initial links (0,1) between nodes given WCM probability
% (2) Assign strengths by distributing all remaining links onto the
% assigned ones in step 1
%
% Returns E, an nxN matrix of all n eigenvalues for all N random modularity
% matrices.
%
% NOTE: if the network has real-valued weights, then the weights are
% converted into multi-edges (i.e. integer weights) by first scaling with
% some factor C; default is 100. This is most appropriate when the weights are e.g.
% correlation values in [0,1]. To omit, or customise, set optional
% conversion factor as outlined below.
%
% ... = POISSONSPARSEWCM(..,C,OPTIONS) sets optional settings:
%       C: sets the conversion factor C; i.e. the amount by which the 
%           weighted adjacency matrix is scaled to get integer weights.
%           C = 'all' sets the conversion factor large enough that the minimum weight
%           is converted to 1. Set C = [] to omit.
%
%       OPTIONS: struct of options:
%           .Expected = {0,1}: if specified (=1), uses the ensemble of generated
%           null models to estimate the expected model. Useful for
%           non-standard null models. [default = 0]
%           .NoLoops = {0,1}: if specified (=1), prevents self-loops in the
%           generated random models [default = 1]
%   
% [..,D,V,X,ALL] = POISSONSPARSEWCM(...) returns:
%           D: a struct, containing diagnostic measurements of the accuracy of the null model 
%           for each of the N repeats, with fields
%           D(i).sAp = strength distribution of the ith repeat
%           D(i).dS = absolute difference between data and ith model strength distributions 
%           D(i).dSN = absolute difference, normalised per node to its strength
%               in the data (i.e. to measure the error relative to magnitude)
%           V: an nxnxN matrix, containing all of the nxn eigenvector matrices of the N repeats            
%           X: the nxn matrix for the expected null model: only
%           returned if Options.Expected = 1;
%           ALL: the nxnxN matrix of every generated network
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
% 13/07/2017: better diagnostics, aligned code with other similiar
% functions
%
% Silvia Maggi & Mark Humphries
addpath('../Helper_Functions/')  % for emptyStruct
addpath('../Network_Analysis_Functions/')  % for expectedA

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

% return all?
blnAll = 0;
if nargout >= 5
    blnAll = 1;
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
P = expectedA(A_int);  

% initialise structs to full storage size, allowing parfor to slice
% appropriately
Pstar = emptyStruct({'Egs','A','V'},[N,1]);

fields = {'SAp','minW','maxW','dS','dSN','Stotal','dStotal','dmax','MDensity','dDensity'};
diagnostics = emptyStruct(fields, [N,1]);

% detect parallel toolbox, and enable if present
blnParallel = autoParallel;

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
    
    ixpairs = find(triu(Aperm,0)>0);     % linear indices of linked pairs for upper triangular matrix
    nLinks = round(Sint/2) - numel(ixpairs);  % number of links left to place: total links - [already placed]
    
    if S ~= K && nLinks > 0 % then is weighted network and there are edges left to place       
        % get as (i,j)
        [irow,jcol] = ind2sub([n,n],ixpairs);
        % get P(link): 
        Plink = sA(irow) .* sA(jcol);
        Plink = Plink ./ sum(Plink); % P(link is placed between each pair)
        
        lambda = nLinks .* Plink; % expected number of links

        Nlink = poissrnd(full(lambda));  % Poisson random number of links made

        Aperm(ixpairs) = Aperm(ixpairs) + Nlink'; % add to existing links
        Aperm = Aperm ./ conversion;  % convert back

        if any(isnan(Aperm(:)))
            keyboard
        end
    end
        
    Aperm = Aperm + Aperm';  % make symmetric
    
    %% diagnostics: how far does random model depart?
    diagnostics(iN) = DiagnosticsOfModelFit(A,Aperm);
    
    %% get eigenvalues 
    if Options.Expected 
        Pstar(iN).A = Aperm;  % store permuted network
    else
        %% go ahead and get eigenvalues
        % P is null model for A, assuming A = P + noise
        % B* = P* - P
        % [Pstar(iN).V,Pstar(iN).Egs] = eig(Aperm - P,'vector');  % not using 'vector' option for backwards compatibility
        [Pstar(iN).V,Egs] = eig(Aperm - P);
        Egs = diag(Egs); % extract vector from diagonal
        [Pstar(iN).Egs,ix] = sort(Egs,'descend'); % ensure eigenvalues are sorted in order
        Pstar(iN).V = Pstar(iN).V(:,ix); % also sort eigenvectors
        
        if blnAll
            Pstar(iN).A = Aperm;  % store permuted network
        end
    end
end

if Options.Expected 
    % generate expected matrix
    Aall = zeros(n);
    for iN = 1:N
      Aall = Aall + Pstar(iN).A;
    end
    ExpWCM = Aall ./ N;
    varargout{3} = ExpWCM;  % return this if asked
    % now compute eigenvalues
    for iN = 1:N
        % [Pstar(iN).V,Pstar(iN).Egs] = eig(Pstar(iN).A - ExpWCM,'vector'); % difference between realisation and Expected model
        [Pstar(iN).V,Egs] = eig(Pstar(iN).A - ExpWCM);
        Egs = diag(Egs); % extract vector from diagonal
        [Pstar(iN).Egs,ix] = sort(Egs,'descend'); % ensure eigenvalues are sorted in order
        Pstar(iN).V = Pstar(iN).V(:,ix); % also sort eigenvectors
    end
end

% now collapse all eigenvalues and vectors into matrix
V = zeros(n,n,N,'single');
E = zeros(n,N);
if blnAll 
    Aall = zeros(n,n,N,'single'); 
else 
    Aall = [];
end
for iN = 1:N
    E(:,iN) = Pstar(iN).Egs;
    V(:,:,iN) = Pstar(iN).V;
    if blnAll Aall(:,:,iN) = Pstar(iN).A; end
end

varargout{1} = diagnostics;
varargout{2} = V;
% 3 is assigned to the expected WCM above
varargout{4} = Aall;