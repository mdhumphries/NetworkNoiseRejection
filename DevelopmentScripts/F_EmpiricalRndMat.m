% clear all 
% close all
% 
% % % **********************
% load Sep2711da01_40_120s_Sxy_Gaussian_1s
% A = Sxyall{1};
% load polblogs
% A = Problem.A;
% N=1;
% nrep = 10;
% % **********************
% % %%
% % % For negative weigthed matrix
% % original = (rand(15)-.5)*3;
% % mi = min(original(:));
% % orig = original - mi;
% % % %%

function [E,varargout] = F_EmpiricalRndMat(A,N)
    
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
    % tic 
    %% Step 2: make weights
    
    xfull = emprand(nonzeros(A(:)),length(nonzeros(Aperm)),1);
    
    ixpairs = find(triu(Aperm,1)>0);
    % get as (i,j)
%     [irow,jcol] = ind2sub([n,n],ixpairs);   
    Aperm(ixpairs) = xfull;

    Aperm = Aperm + Aperm'; %/nrep;
    
    %% check for integer-valued network
    if round(A) == A;
        Aperm = round(Aperm);
    end
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

end

allV = [Pstar.Egs];
allV = allV(:);
varargout{1} = diagnostics;

V = zeros(n,n,N);
E = zeros(n,N);
for iN = 1:N
    E(:,iN) = Pstar(iN).Egs;
    V(:,:,iN) = Pstar(iN).V;
end

varargout{2} = V;

