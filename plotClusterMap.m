function [h,I] = plotClusterMap(W,G,varargin)

% PLOTCLUSTERMAP sorted heatmap of detected clusters
% [H,I] = PLOTCLUSTERMAP(W,G) for a set of N data objects, plots the NxN comparison matrix W
% (e.g. correlation or similarity matrix) as a heatmap grouped into
% clusters C, an N-element array of group membership indexes.
%
% ... = PLOTCLUSTERMAP(...,C,'S') uses the C = nx3 array of colormap values for the
% cluster map. Option 'S' sorts the groups by intrinsic similarity before
% plotting
% 
% Returns:
%    H: the figure handle
%    I: the sorted indices of the nodes in left-to-right (& top-to-bottom)
%
% 28/7/2016: initial version
% 30/8/2016: added colormaps for different situations; plot sorted by
%            similarity
% 23/9/2016: plots within groups in order of similarity; increased resolution of heatmap   
%
% Mark Humphries 28/7/2016

lineColor = [0 0 0];

if nargin > 2 && ~isempty(varargin{1})
    % specify colormap
    cmap = varargin{1}; 
else
    if numel(unique(W)) == 2
        cmap = [1 1 1; 0 0 0];
        lineColor = [0.8 0.3 0.4];
    else
        % make a colormap from C Brewer schemes
        cmap = brewermap(15,'*Blues');  % reversed schemes: light = high
        lineColor = [1 1 1];
    end
end

if nargin > 3
     if varargin{2} == 'S'
        G = [[1:numel(G)]' G];
        [newG,Sgrp,Sin,Sout] = sortbysimilarity(G,W); % sort clusters by similarity (new G)
     %  keyboard  
        % sort by group similarity, and within that by node similarity
        maxG = max(newG(:,2));
        I = []; srt = [];
        for iG = 1:numel(Sgrp)
            % get each group, starting with most similar
            ixNodes = find(newG(:,2) == maxG - (iG-1)); % nodes in this group (original indexes)
            % sort in similarity order
            theseSin = Sin(ixNodes);
            [~,ixSin] = sort(theseSin,'descend');
            % add to group array
            srt = [srt; maxG-(iG-1)*ones(numel(ixNodes),1)];  % sorted descending index of groups
            % store remapped indices
            I = [I; ixNodes(ixSin)]; % original node indices in order
            % keyboard
        end
        
%         % then sort into group order
%         [srt,I] = sort(newG(:,2),'descend'); % so that the most similar group is first

        % keyboard
     end
else
    [srt,I] = sort(G,'descend'); % sort into group index order... (group all 1s, 2s etc)
end

lines = [0; find(abs(diff(srt))==1); numel(srt)]+0.5;  % absolute difference, so it doesn't account

figure
h = imagesc(W(I,I));
colormap(cmap);
axis square
hc = colorbar;

% keyboard
% draw outline box around each cluster
for i=2:numel(lines)
    line([lines(i-1) lines(i) lines(i) lines(i-1); lines(i) lines(i) lines(i-1) lines(i-1)],...
         [lines(i-1) lines(i-1) lines(i) lines(i); lines(i-1) lines(i) lines(i) lines(i-1)],...,
         'Color',lineColor,'LineWidth',1)
end