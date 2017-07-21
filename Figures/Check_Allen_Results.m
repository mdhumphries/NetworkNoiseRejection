%% script to further check Allen data output

clear all; close all;

load('../Networks/Allen_Gene_Leaf.mat');                        % get indices into original data-set
load('../Networks/Allen_mouse_brain_atlas/allen_colours.mat');  % colour-coding of every region

load('../Results/Rejected_Allen_Gene_Leaf.mat');
load('../Results/Clustered_Allen_Gene_Leaf.mat');