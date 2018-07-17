function [E,varargout] = poissonFullWCM(A,N,varargin)

% E = POISSONFULLWCM(A,N) takes the weighted, undirected (n x n) adjacency matrix A and creates N
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
% ... = POISSONFULLWCM(..,C,OPTIONS) sets optional settings:
%       C: sets the conversion factor C; i.e. the amount by which the 
%           weighted adjacency matrix is scaled to get integer weights.
%           C = 'all' sets the conversion factor large enough that the minimum weight
%           is converted to 1. Set C = [] to omit.
%
%       OPTIONS: struct of options:
%           .NoLoops = {0,1}: if specified (=1), prevents self-loops in the
%           generated random models [default = 1]
%   
% [..,D,V,ALL] = POISSONFULLWCM(...) returns:
%           D: a struct, containing diagnostic measurements of the accuracy of the null model 
%           for each of the N repeats, with fields
%           D(i).sAp = strength distribution of the ith repeat
%           D(i).dS = absolute difference between data and ith model strength distributions 
%           D(i).dSN = absolute difference, normalised per node to its strength
%               in the data (i.e. to measure the error relative to magnitude)
%           V: an nxnxN matrix, containing all of the nxn eigenvector matrices of the N repeats            
%           ALL: the nxnxN matrix of every generated network
%
% Notes: 
% (1) assumes A is connected;
% (2) To Do: Add negative weights as separate option: can split into (+) and (-) groups, and assign 
% links, then weights	
%
% ChangeLog:
% 13/07/2018: initial version, cloned from poissonSparseWCM
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
if nargout >= 4
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
nLinks = round(sum(sum(A_int))/2);  % total unique links to place

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
        
    Asamp = zeros(n);
    % self-loops or not
    if Options.NoLoops
        Enode = triu(P,1); % expected number of links between each pair of nodes
    else
        Enode = triu(P,0); % expected number of links between each pair of nodes
    end    
    
    ixpairs = find(Enode>0);                % linear indices of included pairs for upper triangular matrix
    [irow,jcol] = ind2sub([n,n],ixpairs);   % get row, colum version
    
    Plink = sA(irow) .* sA(jcol);       % proportional probability of link between two nodes
    Plink = Plink ./ sum(Plink);    % P(link is placed between each pair): normalised to 1 over all links

    lambda = nLinks .* Plink; % expected number of links between each pair
    
    % keyboard
    
    Nlink = poissrnd(full(lambda));  % Poisson random number of links made
    
    Asamp(ixpairs) = Nlink'; % add to existing links
    Asamp = Asamp ./ conversion;  % convert back

    
%     % keyboard
%     Nlink = poissrnd(full(Enode(Enode > 0)));     % Poisson random number of links made
%     
%     for iL = 1:numel(Nlink)
%         Asamp(irow(iL),jcol(iL)) = Nlink(iL) ./ conversion;  % assign sampled weights, and convert back
%     end
    
    Asamp = Asamp + Asamp';  % make symmetric
    
    %% diagnostics: how far does random model depart?
    diagnostics(iN) = DiagnosticsOfModelFit(A,Asamp);
    
    %% get eigenvalues
    [Pstar(iN).V,Egs] = eig(Asamp - P);  % not using 'vector' option for backwards compatibility
    Egs = diag(Egs); % extract vector from diagonal
    [Pstar(iN).Egs,ix] = sort(Egs,'descend'); % ensure eigenvalues are sorted in order
    Pstar(iN).V = Pstar(iN).V(:,ix); % also sort eigenvectors
    
    % if requesting all networks, then store
    if blnAll
        Pstar(iN).A = A;
    end

end


% now collapse all eigenvalues and vectors into matrix
V = zeros(n,n,N,'single');
E = zeros(n,N);
if blnAll A = zeros(n,n,N,'single'); end

for iN = 1:N
    E(:,iN) = Pstar(iN).Egs;
    V(:,:,iN) = Pstar(iN).V;
    if blnAll A(:,:,iN) = Pstar(iN).A; end
end

varargout{1} = diagnostics;
varargout{2} = V;
varargout{3} = A;