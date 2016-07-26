function [allV,varargout] = WeightedConfigModel(A,N,varargin)

% WEIGHTEDCONFIGMODEL expected eigenvalue distribution for weighted configuration model
% V = WEIGHTEDCONFIGMODEL(A,N) takes the weighted, undirected adjacency matrix A and creates N
% random realisations of the modularity matrix B by randomly generating a 
% null model approximating the configuration model P.
% Returns V, the distribution of all eigenvalues for all random modularity
% matrices.
% 
% V = WEIGHTEDCONFIGMODEL(..,C) sets the conversion factor C; i.e. the amount
% by which the weighted adjacency matrix is scaled to get integer weights.
% C = 'all' sets the conversion factor large enough that the minimum weight
% is converted to 1.
%
% [..,D] = WEIGHTEDCONFIGMODEL(...) returns a struct D, containing diagnostic
% measurements of the accuracy of the null model for each of the N repeats, with fields
%       D(i).sAp = strength distribution of the ith repeat
%       D(i).dS = absolute difference between data and ith model strength distributions 
%       D(i).dSN = absolute difference, normalised per node to its strength
%       in the data (i.e. to measure the error relative to magnitude)
%
% Notes: 
% (1) assumes A is connected;
% (2) Add negative weights as separate option: can split into (+) and (-) groups, and assign % links, then weights	
%
% ChangeLog:
% 17/6/2016: added diagnostics
% 23/6/2016: added conversion scale options; added check for integer
% weights
% 25/7/2016: changed to computation of B* = P* - P as basic model   
%            added Parallel Computing Toolbox support for main loop
%            fixed bug: now returns correct eigenvalues
% 26/7/2016: Full WCM model: two-step generative model	 
%
% Mark Humphries 26/7/2016

n = size(A,1);

sA = sum(A);  % original strength distribution
kA = sum(A>0);  % original degree distribution

% catch errors
sAc = sum(A’); 
if sA ~= sAc
	warning(‘Directed networks will be converted to undirected networks’);
	A = (A + A’) / 2; % convert to undirected
	sA = sum(A); 
	kA = sum(A>0);
end

if any(A < 0)
	error(‘WCM only defined for positive weights - for now…’)
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

% weighted configuration model {expectation]
P = expectedA(A);  CHECK THIS

% initialise structs to full storage size, allowing parfor to slice
% appropriately
Pstar = emptyStruct({'Egs'},[N,1]);

fieldnames = {'conversion','sAp','minW','maxW','dS','dSN','dmax'};
diagnostics = emptyStruct(fieldnames, [N,1]);

% detect parallel toolbox, and enable if present
blnParallel = license('test','Distrib_Computing_Toolbox');

if blnParallel
    nCores = feature('numCores');
    if isempty(gcp('nocreate'))
        parpool('local',nCores-1);  % run on all except 1: don't cripple the machine...
    end
end

parfor iN = 1:N

%% Step 1: create links
K = sum(kA);
m_int = K/2; % so how many loops to match all stubs?

% DO RANDOMISATION BELOW ON ONLY LINKS

%% Step 2: make weights
S = sum(sA);

stubs = sum(A_int);  % how many stubs?
m_int = sum(stubs)/2; % so how many loops to match all stubs?

P_appear = stubs ./ sum(stubs);
C_appear = cumsum(stubs ./ sum(stubs));

% DO RANDOMISATION CONSTRAINED TO LINKED NODE PAIRS

    
    % generate random weighted configuration model
%     Aperm = zeros(n);
%     tic
%     rem_stubs = stubs;
%     for iM = 1:m_int
%         ixs = [0 0];
%         allowed = find(rem_stubs > 0);
%         while ixs(1) == ixs(2)
%             % pick pair at random
%             ixs = ceil(rand(1,2) * numel(allowed));
%         end
%         ix1 = allowed(ixs(1)); ix2 = allowed(ixs(2));
%         rem_stubs(ix1) = rem_stubs(ix1) - 1;
%         rem_stubs(ix2) = rem_stubs(ix2) - 1; 
%         Aperm(ix1,ix2) = Aperm(ix1,ix2) + 1; % edge counting
%         Aperm(ix2,ix1) = Aperm(ix1,ix2);  % symmetry
%         if sum(rem_stubs > 0) == 1    % if only 1 is left, then quit
%             break
%         end
%     end
%     toc
    
    % generate random weighted expected-degree model
    % keyboard
    % slower than a loop...
%     tic
%     X1 = arrayfun(@(x) find(x < C_appear,1),rand(m_int,1)); 
%     toc
    
    % all computation time is taken by this random number generation, due
    % to needing to huge number of random edges and store them - make this faster if possible.
    X1 = discreteinvrnd(P_appear,m_int,1); % source nodes
    
    % X2 = X1(randperm(m_int));  % much faster, but would this strongly
    % bias sampling?? 
    X2 = discreteinvrnd(P_appear,m_int,1); % target nodes

    Atemp = zeros(n);
    for iM = 1:m_int
        Atemp(X1(iM),X2(iM)) = Atemp(X1(iM),X2(iM)) + 1;
    end
    Aperm = Atemp + Atemp'; 

    % convert back...
    Aperm = Aperm ./ conversion;
    
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
    Pstar(iN).Egs = eig(Aperm - P);
    
    % keyboard
end

allV = [Pstar.Egs];
allV = allV(:);
varargout{1} = diagnostics;

