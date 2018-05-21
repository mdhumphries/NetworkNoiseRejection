% script to test Otsu vs k-means on pathological consensus matrices
% matrices from example Ca2+ imaging correlation matrix, from Peron data
% problem appears to be: matrix itself is converging, but the convergence
% detection (k-means) is not working.
%
% So: test Otsu

clear all; close all;

addpath ../Helper_Functions/
addpath ../Network_Analysis_Functions/

load CCons_10.mat
