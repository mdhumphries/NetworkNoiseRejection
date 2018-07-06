function LinkedUnivariateScatterPlots(ax,x,Y,linecolor,varargin)

% LINKEDUNIVERIATESCATTERPLOTS plot univariate scatter of N variables, and link
% LINKEDUNIVERIATESCATTERPLOTS(H,X,Y,L) plots into axis H; given the (m) vector X and (nxm) matrix Y,
% plots each column of Y as a univariate scatter, and links each row of Y
% (assuming they have a common variable). The link are color L (3 element
% color-vector); if L = [], the no links are drawn.
%
% LINKEDUNIVERIATESCATTERPLOTS(...,'Property1',Value1,etc) specifies
% plotting properties, in any order:
%       'MarkerFaceColor': Value should be colormap for different colors
%       for each scatter plot (mx3 array)
%       'MarkerEdgeColor': Value should be colormap for different colors
%       for each scatter plot (mx3 array)
%       'MarkerSize': size of each marker
%       'Linewidth': width of the linking lines
%       'strXlabel': cell array of labels for each scatter plot
% 
% Mark Humphries 13/1/2017

p = inputParser;

[nRows,nCols] = size(Y);

if nCols ~= numel(x)
    error('length of x does not match columns of Y')
end

% defaults
sym = '.-';
if isempty(linecolor) sym = '.'; end

defaultMFColor = ones(nRows,3);
defaultMEColor = zeros(nRows,3);
defaultC = []; defaultE = [];
defaultLWidth = 1;
defaultMsize = 5;
defaultLabels = {};

% parse inputs
if verLessThan('matlab','8.4')  % if before 2014b
    addRequired(p,'ax',@isnumeric);
else
    addRequired(p,'ax',@isobject);
end
addRequired(p,'x',@isnumeric);
addRequired(p,'Y',@isnumeric);
addRequired(p,'linecolor',@isnumeric);
% addOptional(p,'C',defaultC,@isnumeric);
% addOptional(p,'E',defaultE,@isnumeric);
addParameter(p,'MarkerFaceColor',defaultMFColor,@isnumeric)
addParameter(p,'MarkerEdgeColor',defaultMEColor,@isnumeric)
addParameter(p,'MarkerSize',defaultMsize,@isnumeric)
addParameter(p,'Linewidth',defaultLWidth,@isnumeric)
addParameter(p,'strXlabel',defaultLabels,@iscell)

parse(p,ax,x,Y,linecolor,varargin{:});

% finish this if we add error bars later
% LINKEDUNIVERIATESCATTERPLOTS(...,C,E) plots error bars per scatter plot,
% given an (m)-size array of centre points, and (m)x2 array of errors
% [upper; lower];
%
% parse(p,ax,x,Y,linecolor,C,E,varargin{:});
% if numel(p.Results.C) ~= 

% plot
hold on
for iR = 1:nRows
    hp = plot(ax,x,Y(iR,:),sym,'Linewidth',p.Results.Linewidth);  % plot points
    if ~isempty(p.Results.linecolor) 
        set(hp,'Color',p.Results.linecolor);  % plot line color 
    end
    for iC = 1:nCols
        % for each marker, assign colour....
        plot(ax,x(iC),Y(iR,iC),'o','MarkerFaceColor',p.Results.MarkerFaceColor(iC,:),'MarkerEdgeColor',p.Results.MarkerEdgeColor(iC,:),'Markersize',p.Results.MarkerSize)
    end
end

set(ax,'XLim',[x(1)-0.25 x(end)+0.25]);
set(ax,'XTick',x)
set(ax,'XTickLabel',p.Results.strXlabel)

