% clear all; close all
% 
% % analysis parameters
% N = 10;        % repeats of permutation
% % load data
% load Sep2711da01_40_120s_Sxy_Gaussian_1s
% A = Sxyall{1};
% % load('Networks/Lesmis.mat');
% % A = full(Problem.A);


function [E,varargout] = RndPoissonConfigModel(A,N)

% E = RNDPOISSCONFIGMODEL(A,N) takes the weighted, undirected (n x n) adjacency matrix A and creates N
% random realisations of the modularity matrix B by randomly generating a 
% null model approximating the configuration model P.
% Returns E, an nxN matrix of all n eigenvalues for all N random modularity
% matrices.
%
% [..,D,V] = RNDPOISSCONFIGMODEL(...) returns:
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
% (2) To Do: Add negative weights as separate option: can split into (+) and (-) groups, and assign 
% links, then weights	
%
% ChangeLog:
% 31/08/2016: changed from weighted configuration model to random
% poissonian generator model. Modified from WEIGHTEDCONFIGMODEL, Mark 
% Humphries's script 26/07/2016
%
% Silvia Maggi 31/08/2016

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

% weighted configuration model [expectation]
P = expectedA(A);  

% initialise structs to full storage size, allowing parfor to slice
% appropriately
Pstar = emptyStruct({'Egs','A','V'},[N,1]);

fieldnames = {'sAp','minW','maxW','dS','dSN','dmax','Aperm'};
diagnostics = emptyStruct(fieldnames, [N,1]);

% detect parallel toolbox, and enable if present
blnParallel = license('test','Distrib_Computing_Toolbox');

if blnParallel
    nCores = feature('numCores');
    if isempty(gcp('nocreate'))
        parpool('local',nCores);  % run on all
    end
end


parfor iN = 1:N;
% for iN = 1:N
    %% Step 1: create links
    K = sum(kA);  % total number of links
   
    pnode = expectedA(A>0); % probability of link between each pair of nodes
    
    Aperm = real(rand(n,n) < triu(pnode,1));  % add link to all that pass test
    
    
    %% Step 2: make weights
    S = sum(sA);  % total weights
    
%     if S ~= K  % then is weighted network  
        rndmat = zeros(size(A));
        % disp('Weighted network')
        % linear indices of linked pairs for upper triangular matrix
        ixpairs = find(triu(Aperm,1)>0); 
        
        % get as (i,j)
        [irow,jcol] = ind2sub([n,n],ixpairs);
        % linear indices of linked pairs for lower triangular matrix
        % (symmetric respect to the diagonal)
        linearInd = sub2ind([n,n], jcol, irow);
        % get P(link) in order.... 
        Plink = poissrnd(sA(irow) .* sA(jcol));
        
        rndmat(ixpairs) = Plink;
        rndmat(linearInd) = Plink;

        kkk = rndmat.*(ones(size(A,1))-eye(size(A,1))); % remove diagonal
        Aperm = sum(sum(A))/sum(sum(kkk))*kkk; % normalisation
        

    %% diagnostics: how far does random model depart?
    diagnostics(iN).sAp = sum(Aperm);  % degree
    diagnostics(iN).minW = min(Aperm); % minimum weight
    diagnostics(iN).maxW = max(Aperm);    % maximum weight
    
    % figure; ecdf(sA); hold on; ecdf(sAp); title('Degree distributions of original and permuted network')
    
    diagnostics(iN).dS = abs(sA - diagnostics(iN).sAp);
    diagnostics(iN).dSN = 100* diagnostics(iN).dS ./ sA; % difference as fraction of original degree
 
    diagnostics(iN).dmax =  max(A) - diagnostics(iN).maxW;
    diagnostics(iN).Aperm = Aperm;
    % figure; ecdf(dSN); title('ECDF of error as proportion of original degree')
    
    %% get eigenvalues
    % P is null model for A, assuming A = P + noise
    % B* = P* - P
    [Pstar(iN).V,Pstar(iN).Egs] = eig(Aperm - P,'vector');
    [Pstar(iN).Egs,ix] = sort(Pstar(iN).Egs,'descend'); % ensure eigenvalues are sorted in order
    Pstar(iN).V = Pstar(iN).V(:,ix); % also sort eigenvectors
%     Pstar(iN).A = Aperm;
    % keyboard
end

varargout{1} = diagnostics;

V = zeros(n,n,N);
E = zeros(n,N);
for iN = 1:N
    E(:,iN) = Pstar(iN).Egs;
    V(:,:,iN) = Pstar(iN).V;
end

varargout{2} = V;