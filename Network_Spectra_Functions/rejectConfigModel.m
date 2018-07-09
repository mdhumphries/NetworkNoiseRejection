function [Control,varargout] = rejectConfigModel(A,pars,varargin)

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
% [C,R] = REJECTCONFIGMODEL(...,'Reject',RPARS) also runs the node
% rejection on the low-D projection, with parameters RPARS, returning
% results in struct R, and adding the following to C:
%       .Asignal_final : the final weight matrix, with nodes rejected and
%                               leaves stripped
%       .ixSignal_Final : IDs of nodes retained in the final weight matrix  
%
%  09/07/2018 : initial version
%
% Mark Humphries 

addpath('../Network_Analysis_Functions/')  % for prep_A

blnReject = 0;

if nargin > 2
    if varargin{1} == 'Reject'; blnReject = 1; end
    optionsReject = varargin{2};
end
    
% Prepare A for Noise Rejection
[Data.A,Data.ixRetain,Data.Comps,Data.CompSizes] = prep_A(A);

Control.P = expectedA(Data.A);  % standard configuration model expectation
B = Data.A - Control.P;         % modularity matrix

% compute based on spectral rejection
[Control.Emodel,~,Vmodel] = RndPoissonConfigModel(Data.A,pars.N,pars.C);  % sample from space of configuration models
[Control.Dspace,~,Control.Dn,Control.EigEst] = LowDSpace(B,Control.Emodel,pars.alpha); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors

% compute groups based on just positive eigenvalues
egs = eig(B);  % eigenspectra of data modularity matrix
egs = sort(egs,'descend'); % sort eigenvalues into descending order 
Control.PosDn = sum(egs > pars.eg_min);

% reject nodes if also requested
if blnReject
    Rejection = NodeRejection(B,Control.Emodel,pars.alpha,Vmodel,optionsReject); % N.B. also calls LowDSpace function to find projections
        varargout{1} = Rejection;

        
    if Control.Dn > 0 
        % new signal matrix
        Control.Asignal = Data.A(Rejection.ixSignal,Rejection.ixSignal);

        % connected signal matrix: find largest component, and use that - store
        % others
        [Control.Asignal_comp,ixRetain,Control.SignalComps,Control.SignalComp_sizes] = prep_A(Control.Asignal); 
        Control.ixSignal_comp = Rejection.ixSignal(ixRetain);  % original node indices

        % and then strip out leaves - nodes with single links
        K = sum(Control.Asignal_comp);
        ixLeaves = find(K==1); ixKeep = find(K > 1);

        Control.ixSignal_Final = Control.ixSignal_comp(ixKeep);
        Control.ixSignal_Leaves = Control.ixSignal_comp(ixLeaves);
        Control.Asignal_final = Control.Asignal_comp(ixKeep,ixKeep);
    else
        Control.Asignal = []; Control.Asignal_final = []; Control.Asignal_comp = [];
        Control.ixSignal_Final = [];
        Control.ixSignal_Leaves = [];

    end
end