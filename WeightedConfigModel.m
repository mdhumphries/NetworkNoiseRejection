function [E,varargout] = WeightedConfigModel(A,N,varargin)

% WEIGHTEDCONFIGMODEL expected eigenvalue distribution for weighted configuration model
% E = WEIGHTEDCONFIGMODEL(A,N) takes the weighted, undirected (n x n) adjacency matrix A and creates N
% random realisations of the modularity matrix B by randomly generating a 
% null model approximating the configuration model P.
% Returns E, an nxN matrix of all n eigenvalues for all N random modularity
% matrices.
% 
% ... = WEIGHTEDCONFIGMODEL(..,C) sets the conversion factor C; i.e. the amount
% by which the weighted adjacency matrix is scaled to get integer weights.
% C = 'all' sets the conversion factor large enough that the minimum weight
% is converted to 1.
%
% [..,D,V] = WEIGHTEDCONFIGMODEL(...) returns:
%           a struct D, containing diagnostic measurements of the accuracy of the null model 
%           for each of the N repeats, with fields
%           D(i).sAp = strength distribution of the ith repeat
%           D(i).dS = absolute difference between data and ith model strength distributions 
%           D(i).dSN = absolute difference, normalised per node to its strength
%               in the data (i.e. to measure the error relative to magnitude)
%           an nxnxN matrix V, containing all of the nxn eigenvector matrices of the N repeats            
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
%
% Mark Humphries 26/7/2016

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

% values for diagnostic checking
minW = min(min(A(A>0)));
maxW = max(max(A));

% quantisation steps
if nargin >= 3
    conversion = varargin{1};
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

fieldnames = {'conversion','sAp','minW','maxW','dS','dSN','dmax'};
diagnostics = emptyStruct(fieldnames, [N,1]);

% detect parallel toolbox, and enable if present
% blnParallel = license('test','Distrib_Computing_Toolbox');
% 
% if blnParallel
%     nCores = feature('numCores');
%     if isempty(gcp('nocreate'))
%         parpool('local',nCores);  % run on all
%     end
% end


% parfor iN = 1:N
for iN = 1:N
    %% Step 1: create links
    K = sum(kA);  % total number of links
   
    pnode = expectedA(A>0); % probability of link between each pair of nodes
    
    Aperm= real(rand(n,n) < triu(pnode,1));  % add link to all that pass test
    
    
    %% Step 2: make weights
    S = sum(sA);  % total weights
    
    sAint = sum(A_int); % integer strength
    Sint = sum(sAint); % integer total strength
    if S ~= K  % then is weighted network    
        
        ixpairs = find(triu(Aperm,1)>0);   % linear indices of linked pairs
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
    diagnostics(iN).conversion = conversion; % store
    diagnostics(iN).sAp = sum(Aperm);  % degree
    diagnostics(iN).minW = min(Aperm); % minimum weight
    diagnostics(iN).maxW = max(Aperm);    % maximum weight
    
 
    % figure; ecdf(sA); hold on; ecdf(sAp); title('Degree distributions of original and permuted network')
    
    diagnostics(iN).dS = abs(sA - diagnostics(iN).sAp);
    diagnostics(iN).dSN = 100* diagnostics(iN).dS ./ sA; % difference as fraction of original degree
 
    diagnostics(iN).dmax =  max(A) - diagnostics(iN).maxW;
    % figure; ecdf(dSN); title('ECDF of error as proportion of original degree')
    
    %% get eigenvalues
    % P is null model for A, assuming A = P + noise
    % B* = P* - P
    [Pstar(iN).V,Pstar(iN).Egs] = eig(Aperm - P,'vector');
    % Pstar(iN).A = Aperm;
    % keyboard
end

% Aall = zeros(n);
% for iN = 1:N
%   Aall = Aall + Pstar(iN).A;
% end
% ExpWCM = Aall ./ N;
% keyboard

% E = [Pstar.Egs];
% E = E(:);
varargout{1} = diagnostics;

V = zeros(n,n,N);
E = zeros(n,N);
for iN = 1:N
    E(:,iN) = Pstar(iN).Egs;
    V(:,:,iN) = Pstar(iN).V;
end

varargout{2} = V;

