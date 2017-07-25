% Allen_colours_imageplot.m
% 
% Script to make an image of the sorted-by-cluster 'allen_colours'
%
% M. Evans July 2017

clear all; close all

% Change path as appropriate
load('/Users/mathew/work/NetworkNoiseRejection/Networks/Allen_Gene_Leaf.mat')
load('/Users/mathew/work/NetworkNoiseRejection/Networks/Allen_mouse_brain_atlas/allen_colours.mat')
load('/Users/mathew/work/NetworkNoiseRejection/Results/Clustered_Allen_Gene_Leaf.mat', 'Full')

% Qmax
figure;
[~,ix] = sort(Full.QmaxCluster);
imagesc(repmat(1:numel(ix),6,1))
hold all
plot(Full.QmaxCluster(ix),'k','linewidth',2)

colormap(double(allen_colours(leafIDs_final(ix),:))/255)

title('Qmax')
ylabel('Cluster assignment')
xlabel('Sorted region colours')
axis square

% Consensus clustering
figure;
[~,ix2] = sort(Full.ConsCluster);
imagesc(repmat(1:numel(ix),27,1))
hold all
plot(Full.ConsCluster(ix2),'k','linewidth',2)

colormap(double(allen_colours(leafIDs_final(ix2),:))/255)

title('Consensus')
ylabel('Cluster assignment')
xlabel('Sorted region colours')
axis square
