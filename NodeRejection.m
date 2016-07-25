function [D,varargout] = NodeRejection(B,allV,bnds,varargin)

% NODEREJECTION separates nodes into "signal" and "noise"
% D = NODEREJECTION(B,V,I) splits the nodes in a network into signal 
% and noise components, given: 
%       B: the modularity matrix of the network, defined using a null model (e.g Weighted Configuration Model)
%       V: the null-model eigenvalue distribution (from e.g. expectedEigsUnd) 
%       I: specified rejection interval (propotion: 0.05, for 95%; 0.01 for
%       99%, and so on); if I is specified as an n-length array {I1,I2,...,In], 
%       then a decompositin will be returned for each I  
%
% Returns: D, an n-length struct array with fields:
%               .ixSignal: the node indices in the signal component of the
%               network
%               .ixNoise: the node indices in the noise component of the
%               network
%                           
% Notes: 
% (1) assumes A is connected;
%
%
% ChangeLog:
% 25/7/2016: initial version
%
% Mark Humphries 25/7/2016

% compute eigenvalues & vectors of modularity matrix
[V,D] = eig(B);
egs = diag(D);

% find bounds, and calculate numbers to reject

% project and reject

