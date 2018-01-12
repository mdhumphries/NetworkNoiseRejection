function exportPPTfig(h,figID,varargin)

% EXPORTPPTFIG quick export of raw MATLAB figures for presentation
% EXPORTPPTFIG(H,N) exports the figure with handle H to the filename in
% string N. Saves the path hard-coded in this function
% 
% EXPORTPPTFIG(...,S) sets the figure to size [L,B,W,H]: left L and
% bottom-up B on screen, and W wide and H high (in cm)
%
% exportpath = 'C:\Users\mqbssmhg.DS\Dropbox\My Talks\Manchester\Data Clubs\Mar16_2017\';
exportpath = '/Users/mqbssmhg/Dropbox/My Papers/Conference Abstracts/NetSci2018/';

% common parameters
format = 'png'; %'tiffn' for submission
color = 'rgb';
dpi = 800;
fontsize = 7;
linewidth = 1;
boxlinewidth = 0.5;
M = 2;  % marker size
flag = '3'; % clusters only

if nargin > 2
    PPTsize = varargin{1};
else
    PPTsize = [5 5 6 6];
end


orig_Xlim = get(gca,'XLim');
orig_Ylim = get(gca,'YLim');

set(h,'Units', 'centimeters', 'PaperPositionMode', 'auto','Position',PPTsize); 
% keyboard

set(gca,'XLim',orig_Xlim);
set(gca,'YLim',orig_Ylim);

% keyboard

htext = findobj(gcf, 'type', 'text');
set(htext,'FontName','Helvetica','FontSize',fontsize);

hAxes = findobj(gcf, 'type', 'axes');
for i=1:numel(hAxes)
    set(get(hAxes(i),'XLabel'),'FontName','Helvetica','FontSize',fontsize);
    set(get(hAxes(i),'YLabel'),'FontName','Helvetica','FontSize',fontsize);
    set(hAxes(i),'FontName','Helvetica','FontSize',fontsize);
    set(hAxes(i),'Box','off','TickDir','out','LineWidth',boxlinewidth);
end


exportfig(gcf,[exportpath 'Fig_' figID],'Color',color,'Format',format,'Resolution',dpi)
