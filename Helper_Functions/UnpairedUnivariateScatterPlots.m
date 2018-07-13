function UnpairedUnivariateScatterPlots(ax,Y,varargin)

% UNPAIREDUNIVERIATESCATTERPLOTS plot univariate scatter of N unpaired variables
% UNPAIREDUNIVERIATESCATTERPLOTS(H,Y) plots into axis H; given the m-length cell array Y,
% plots each cell of Y as a univariate scatter, in order of the array
%
% UNPAIREDUNIVERIATESCATTERPLOTS(...,C,E) plots error bars per scatter plot,
% given an (m)-size array of centre points C, and (m)x2 array of errors E =
% [upper; lower]; plots an horizontal line at C, and a line between C+upper
% and C-lower.
%
% UNPAIREDUNIVERIATESCATTERPLOTS(...,'Property1',Value1,etc) specifies
% plotting properties, in any order:
%       'MarkerFaceColor': Value should be colormap for different colors
%       for each scatter plot (mx3 array)
%       'MarkerEdgeColor': Value should be colormap for different colors
%       for each scatter plot (mx3 array)
%       'strXlabel': cell array of labels for each scatter plot
%       'ErrorBarColor': Value should be colormap for different colors
%       for each scatter plot (mx3 array)
%       'ErrorBarWidth': line width (in px)
%       'ErrorBarPosition': {'front'} | 'back': defines where to plot the error-bars with respect to the scatters.   
%
% Mark Humphries 16/11/2016

p = inputParser;

[nData] = numel(Y);
x = 1:nData;

% defaults
sym = '.';

defaultMFColor = ones(nData,3);
defaultMEColor = zeros(nData,3);
defaultC = []; defaultE = [];
defaultMsize = 5;
defaultLabels = {};
defaultPosition = 'front';
defaultEWidth = 1;
defaultErrColor = repmat([0.6 0.6 0.6],nData,1);

% parse inputs
if verLessThan('matlab','8.4')  % if before 2014b
    addRequired(p,'ax',@isnumeric);
else
    addRequired(p,'ax',@isobject);
end
addRequired(p,'Y',@iscell);
addOptional(p,'C',defaultC,@isnumeric);
addOptional(p,'E',defaultE,@isnumeric);
addParameter(p,'MarkerFaceColor',defaultMFColor,@isnumeric);
addParameter(p,'MarkerEdgeColor',defaultMEColor,@isnumeric);
addParameter(p,'MarkerSize',defaultMsize,@isnumeric);
addParameter(p,'strXlabel',defaultLabels,@iscell);
addParameter(p,'ErrorBarPosition',defaultPosition,@ischar);
addParameter(p,'ErrorBarColor',defaultErrColor,@isnumeric);
addParameter(p,'ErrorBarWidth',defaultEWidth,@isnumeric);


% parse(p,ax,Y,varargin{:});
parse(p,ax,Y,varargin{:});

if ~isempty(p.Results.C) && numel(p.Results.C) ~= numel(p.Results.Y)
    error('Number of data-sets and centre-points do not agree');
end

if ~isempty(p.Results.C) && numel(p.Results.C) ~= numel(p.Results.Y)
    error('Number of data-sets and centre-points do not agree');
end

if ~isempty(p.Results.E)
    E = p.Results.E;
    if size(p.Results.E,2) > 2   % re-orient to get 2 columns
        E = E';
    end
end

% plot
hold on

for iC = 1:nData
    if strcmp(p.Results.ErrorBarPosition,'front')
        plot(ax,x(iC),Y{iC},'o','MarkerFaceColor',p.Results.MarkerFaceColor(iC,:),'MarkerEdgeColor',p.Results.MarkerEdgeColor(iC,:),'Markersize',p.Results.MarkerSize)
    end
    if ~isempty(p.Results.C)
        % error bars
        line([x(iC)-0.25 x(iC)+0.25],[p.Results.C(iC), p.Results.C(iC)],'Color',p.Results.ErrorBarColor(iC,:),'Linewidth',p.Results.ErrorBarWidth);
        line([x(iC) x(iC)],[p.Results.C(iC)+E(iC,1),p.Results.C(iC)-E(iC,2)],'Color',p.Results.ErrorBarColor(iC,:),'Linewidth',p.Results.ErrorBarWidth);
    end
    
    if strcmp(p.Results.ErrorBarPosition,'back')
        plot(ax,x(iC),Y{iC},'o','MarkerFaceColor',p.Results.MarkerFaceColor(iC,:),'MarkerEdgeColor',p.Results.MarkerEdgeColor(iC,:),'Markersize',p.Results.MarkerSize)
    end     
end


set(ax,'XLim',[x(1)-0.25 x(end)+0.25]);
set(ax,'XTick',x)
set(ax,'XTickLabel',p.Results.strXlabel)

