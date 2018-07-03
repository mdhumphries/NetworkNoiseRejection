function drawgrid(h,x,y,xlim,ylim,colour,width)

% DRAWGRID draws a line grid on 2D axes, e.g. imagesc
% DRAWGRID(H,X,Y,XL,YL,C,W) draws a grid on axes H, defined by:
%       X : an array of X values for vertical lines
%       Y : an array of Y values for horiztonal lines
%       XL: the [L,U] array of (L)ower and (U)pper limits of horizontal
%       lines
%       YL: the [L,U] array of (L)ower and (U)pper limits of vertical
%       lines
%       C : a 3 element array of line colour
%       W : a scalar for line width (pts)
%
% 03/07/2018: initial version
%
% Mark Humphries

axes(h); % select specified axes
for i = 1:numel(x)
    line([x(i) x(i)],ylim,'Color',colour,'Linewidth',width);
end

for i = 1:numel(y)
    line(xlim,[y(i) y(i)],'Color',colour,'Linewidth',width);
end 