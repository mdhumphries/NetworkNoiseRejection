function [E,varargout] = linkSparseWCM(A,N,varargin)

% LINKSPARSEWCM expected eigenvalue distribution for sparse weighted configuration model
% E =  LINKSPARSEWCM(A,N) takes the weighted, undirected (n x n) adjacency matrix A and creates N
% random realisations of the modularity matrix B by randomly generating a 
% null model approximating the configuration model P.
% Returns E, an nxN matrix of all n eigenvalues for all N random modularity
% matrices.
% 
% NOTE: if the network has real-valued weights, then the weights are
% converted into multi-edges (i.e. integer weights) by first scaling with
% some factor C; default is 100. This is most appropriate when the weights are e.g.
% correlation values in [0,1]. To omit, or customise, set optional
% conversion factor as outlined below.
%
% ... =  LINKSPARSEWCM(..,C,OPTIONS) sets optional settings:
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
% [..,D,V,X,ALL] = LINKSPARSEWCM(...) returns:
%           D: a struct, containing diagnostic measurements of the accuracy of the null model 
%           for each of the N repeats, with fields
%           D(i).sAp = strength distribution of the ith repeat
%           D(i).dS = absolute difference between data and ith model strength distributions 
%           D(i).dSN = absolute difference, normalised per node to its strength
%               in the data (i.e. to measure the error relative to magnitude)
%           V: an nxnxN matrix, containing all of the nxn eigenvector matrices of the N repeats            
%           X: the nxn matrix for the expected configuration model: only
%           returned if Options.Expected = 1;%
%           ALL: the nxnxN matrix of every generated network
%
% Notes: 
% (1) assumes A is connected;
% (2) To Do: Add negative weights as separate option: can split into (+) and (-) groups, and assign % links, then weights	
%
% ChangeLog:
% 17/6/2016: added diagnostics
% 23/6/2016: added conversion scale options; added check for integer
% weights
% 25/7/2016: changed to computation of B* = P* - P as basic model   
%            added Parallel Computing Toolbox support for main loop
%            fixed bug: now returns correct eigenvalues
% 26/7/2016: Full WCM model: two-step generative model	 
%            Returns eigenvectors for all generated null models
% 30/8/2016: Option to return expected network
%            Option to eliminate self-loops
%
% Mark Humphries 

addpath('../Helper_Functions/')  % for emptyStruct and discreteinvrnd
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

% values for diagnostic checking
minW = min(min(A(A>0)));
maxW = max(max(A));

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
   
    % use adjacency matrix form: (A>0); define probability of a link based
    % on node degree (standard configuration model)
    pnode = expectedA(A>0); % probability of link between each pair of nodes, defined only on existing links
     
    % test those probabilities to create links; with option to ignore self-loops
    if Options.NoLoops
        Aperm = real(rand(n,n) < triu(pnode,1));  % don't include diagonal: no self-loops (slightly underestimates degree)
    else
        Aperm = real(rand(n,n) < triu(pnode,0));  % include diagonal: allow self-loops
    end
    
    % keyboard
    %% Step 2: make weights
    S = sum(sA);  % total weights
    
    sAint = sum(A_int); % integer strength
    Sint = sum(sAint); % integer total strength
    if S ~= K  % then is weighted network    
        
        % disp('Weighted network')
        ixpairs = find(triu(Aperm,0)>0);   % linear indices of linked pairs
        % get as (i,j)
        [i,j] = ind2sub([n,n],ixpairs);
        
        % get P(link) in order.... 
        Plink = sA(i) .* sA(j);
        Plink = Plink ./ sum(Plink); % P(link is placed between each pair)
        
        % randomly generate pairs of links
        nLinks = round(Sint/2) - numel(ixpairs);  % total links - [already placed]
        
        % keyboard
        
        X1 = discreteinvrnd(Plink,nLinks,1); % randomly sampled indices of pairs
        
        % add up the links
        for iM = 1:nLinks
            Aperm(ixpairs(X1(iM))) = Aperm(ixpairs(X1(iM))) + 1;
        end
             
        % convert back...
        Aperm = Aperm ./ conversion;
    end
    
    % complete network: undirected
    Aperm = Aperm + Aperm'; 

    %% diagnostics: how far does random model depart?
    diagnostics(iN) = DiagnosticsOfModelFit(A,Aperm);
    
    if Options.Expected 
        Pstar(iN).A = Aperm;  % store permuted network
    else
        %% go ahead and get eigenvalues
        % P is null model for A, assuming A = P + noise
        % B* = P* - P
%         [Pstar(iN).V,Pstar(iN).Egs] = eig(Aperm - P,'vector');
%         [Pstar(iN).Egs,ix] = sort(Pstar(iN).Egs,'descend'); % ensure eigenvalues are sorted in order
        [Pstar(iN).V,Egs] = eig(Aperm - P);  % not using 'vector' option for backwards compatibility
        Egs = diag(Egs); % extract vector from diagonal
        [Pstar(iN).Egs,ix] = sort(Egs,'descend'); % ensure eigenvalues are sorted in order
        Pstar(iN).V = Pstar(iN).V(:,ix); % also sort eigenvectors
        if blnAll
            Pstar(iN).A = Aperm;  % store permuted network
        end
    end
   
    % keyboard
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
%         [Pstar(iN).V,Pstar(iN).Egs] = eig(Pstar(iN).A - ExpWCM,'vector');
%         [Pstar(iN).Egs,ix] = sort(Pstar(iN).Egs,'descend'); % ensure eigenvalues are sorted in order
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
