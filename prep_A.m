% Function to go through all the steps to prepare a connectivity matrix for
% noise rejection. Namely, 

function [newA,nz_e] = prep_A(A)
% Ensure A is dense, not sparse
A = full(A);

% Keep only non-zero elements
nz_e = find(sum(A)); % nonzero_elements

A = A(nz_e,nz_e);

% Remove diagonal elements
A(find(eye(length(A)))) = 0;

newA = A;


