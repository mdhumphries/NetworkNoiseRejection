% High level function to run noise rejection on an example network, and
% save out an array of signal (1) and noise (0) labels. 
%
% Input: A - connectivity matrix
%        filename - filename to use for saving signal node array

function signal_nodes = reject_the_noise(A,filename)

% Prepare A for Noise Rejection
[newA,nz_e] = prep_A(A);

N = ceil(20000/length(A));        % repeats of permutation. Aiming to get ~20000 samples.
alpha = 0;      % confidence interval on estimate of maxiumum eigenvalue for null model

% WCM model options
WCMOptions.Expected = 1;
WCMOptions.NoLoops = 1;

% NodeRejection options
options.Weight = 'linear'; % 'linear' is default
options.Norm = 'L2'; % L2 is default

% get expected distribution of eigenvalues under null model (here, WCM)
[Emodel,diagnostics,Vmodel,ExpWCM] = WeightedConfigModel(newA,N,'all',WCMOptions);

% Decompose into signal and noise
B = newA - ExpWCM;  % modularity matrix using chosen null model

% compare data and model
Edata = eig(B);

% find low-dimensional projection
[Dspace,Ix,Dn,EigEst] = LowDSpace(B,Emodel,alpha); % to just obtain low-dimensional projection

% node rejection within low-dimensional projection
R = NodeRejection(B,Emodel,alpha,Vmodel,options); % N.B. also calls function to find projections

%% Print out array of signal(1) and noise(0) nodes to .csv

signal_nodes = zeros(length(A),1);
signal_nodes(nz_e(R.ixSignal)) = 1;

save(['Networks/signal_',filename,'.mat'],'signal_nodes')

