function h = plotClusterMap(W,G,varargin)

% PLOTCLUSTERMAP sorted heatmap of detected clusters
% H = PLOTCLUSTERMAP(W,G) for a set of N data objects, plots the NxN comparison matrix W
% (e.g. correlation or similarity matrix) as a heatmap grouped into
% clusters C, an N-element array of group membership indexes.
%
% 28/7/2016: initial version
%
% Mark Humphries 28/7/2016

if nargin > 2
    % specify colormap
    colormap(varargin{1});
end

[srt,I] = sort(G,'ascend');
lines = [0; find(diff(srt)==1); numel(G)]+0.5;

figure
h = imagesc(W(I,I));

% draw outline box around each cluster
for i=2:numel(lines)
    line([lines(i-1) lines(i) lines(i) lines(i-1); lines(i) lines(i) lines(i-1) lines(i-1)],...
         [lines(i-1) lines(i-1) lines(i) lines(i); lines(i-1) lines(i) lines(i) lines(i-1)],...,
         'Color',[1 1 1],'LineWidth',1)
end