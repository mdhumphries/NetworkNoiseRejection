function Control = rejectConfigModel(A,pars)

% REJECTCONFIGMODEL do spectral rejection on original configuration model
% C = REJECTCONFIGMODEL(W,PARS) does spectral rejection on the undirected,
% weighted network W using the original configuration model as the null
% model. Struct PARS has fields:
%       .N :    repeats of permutation
%       .alpha : confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
%       .C   : conversion factor for real-valued weights (set=1 for integers)
%       .eg_min : given machine error, what is acceptable as "zero" eigenvalue
%
% Returns struct C, with fields:
%       .P : the expectation of the null model for the data network
%       .Dn : the number N of retained dimensions
%       .Dspace : the retained low-dimensional space (top N eigenvectors)
%       .PosDn: the number of retained dimensions is just using positive
%       eigenvalues
%
%  09/07/2018 : initial version
%
% Mark Humphries 

addpath('../Network_Analysis_Functions/')  % for prep_A

% Prepare A for Noise Rejection
[Data.A,Data.ixRetain,Data.Comps,Data.CompSizes] = prep_A(A);

Control.P = expectedA(Data.A);  % standard configuration model expectation
B = Data.A - Control.P;         % modularity matrix

% compute based on spectral rejection
[Control.Emodel,~,~] = RndPoissonConfigModel(Data.A,pars.N,pars.C);  % sample from space of configuration models
[Control.Dspace,~,Control.Dn,Control.EigEst] = LowDSpace(B,Control.Emodel,pars.alpha); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors

% compute groups based on just positive eigenvalues
egs = eig(B);  % eigenspectra of data modularity matrix
egs = sort(egs,'descend'); % sort eigenvalues into descending order 
Control.PosDn = sum(egs > pars.eg_min);
