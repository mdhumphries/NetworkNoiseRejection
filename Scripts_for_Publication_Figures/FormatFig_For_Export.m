function FormatFig_For_Export(h,fontsize,fontname,boxlinewidth)

% FORMATANDEXPORTFIG does what it says on the tin
% FORMATANDEXPORTFIG(H,FONTSIZE,FONTNAME,BOXLINEWIDTH)
% Given figure handle H, formats the figure's text, ticks and boxes to
% defaults for exporting
%
% MDH 8/2/2018

%% formatting
htext = findobj(h, 'type', 'text');
set(htext,'FontName',fontname,'FontSize',fontsize);

hAxes = findobj(h, 'type', 'axes');
for i=1:numel(hAxes)
    set(get(hAxes(i),'XLabel'),'FontName',fontname,'FontSize',fontsize);
    set(get(hAxes(i),'YLabel'),'FontName',fontname,'FontSize',fontsize);
    set(hAxes(i),'FontName',fontname,'FontSize',fontsize);
    set(hAxes(i),'Box','off','TickDir','out','LineWidth',boxlinewidth);
end
