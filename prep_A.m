% Function to go through all the steps to prepare a connectivity matrix for
% noise rejection. Namely, 

function [A,nz_e] = prep_A(A)

% Keep only non-zero elements
nz_e = find(sum(A)); % nonzero_elements

A = A(nz_e,nz_e);

% Remove diagonal elements
A((eye(length(A)))==1) = 0;



