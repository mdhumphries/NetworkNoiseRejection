function [h,xpad,ypad,xtick,ytick] = plotMatrix(X,Y,M,L,varargin)

% PLOTMATRIX plots a checkerboard plot of a matrix
% [H,Xplot,Yplot,Xtick,Ytick] = PLOTMATRIX(X,Y,M,L) plots matrix M as a checkerboard plot, using all
% rows and columns, aligned to axes X and Y: 
%   M : n x m matrix of scalars
%   X : n-length array of values
%   Y : m-length array of values
%   L : [x y] the lower limit of the X and Y scales (e.g. 0); all rows and
%   columns are sized by the difference between consecutive values in X 
%   (columns) and Y (rows); the lower limit sets the width of the first
%   row and column
%   
% ... = PLOTMATRIX(...,C): optional argument specifying colourmap
%
% Returns:
%   H : handle to the PCOLOR plot
%   Xplot, Yplot: the arrays used to define the plotted x and y axes
%   (useful for e.g. using DRAW_GRID)
%   Xtick, Ytick: the arrays used to define the location of the axis tick
%   marks
%
% NOTES:
%   Essentially this solves the problem that PCOLOR cuts off the last row
%   and column; and that IMAGE (and IMAGESC) cannot handle unequal sized
%   axes. This function uses PCOLOR, with adjusted matrix and axes to plot
%   the full matrix
% 
% 05/07/2018: initial version
%
% Mark Humphries

cmap = [];
if nargin > 4
    cmap = varargin{1};
end
    
plotM = padarray(M,[1 1],'replicate','post');
xpad = [L(1) X];
ypad = [L(2) Y];
xtick = xpad(1:end-1) + diff(xpad)/2;
ytick = ypad(1:end-1) + diff(ypad)/2;

h = pcolor(xpad,ypad,plotM);
set(gca,'XTick',xtick,'XTickLabel',X,'YTick',ytick,'YTickLabel',Y)

if cmap
    colormap(cmap);
end
