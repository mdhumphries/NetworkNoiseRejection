function graphplot2d(xy,W,varargin)
% This program was written to take in xy-coordinates of a network as well as
% the adjacency matrix W, to graph the network using these coordinates.  We
% have many programs that can calculate these coordinates either respecting
% community structure (KamadaC.m, fruc_reinC.m, KKFR.m, FRKK.m) or ignoring 
% it (Kamada.m, fruc_rein.m).  This program can also take in a factor for edge length, 
% colors for each node and shapes for each node.  If colors and alpha are given only,
% this program sets up the skeleton for a very nice legend.
%
%Inputs:
%
% xy= matrix of xy-coordinates for each node in the network
%
% W=Adjacency Matrix for the network
%
% Optional Variables:
% 
% alpha = factor for edge colouring
%
% colors = a vector of numbers defining colors, if zero is one of the
% colors, that is given as purple; or N-by-3 matrix specifying a color [R G B] for
% each node
%
% shapes = a vector of strings defining the shapes for each node, i. e. 'd'
% for diamond, '.' for dot, etc.; set
%
% size = a scalar setting the size of the node shape(s)
%
% Function Calls:
%
% graphplod2d(xy,W)
% graphplot2d(xy,W,alpha)
% graphplot2d(xy,W,alpha,colors)
% graphplot2d(xy,W,alpha,colors,shapes)
% graphplot2d(xy,W,alpha,colors,shapes,size)
%
% Last Modified by ALT May 19,2009, Created by PJM


if (length(varargin)) & ~isempty(varargin{1})
    alpha=varargin{1};
else
    alpha=2;
end

%
%Set colormap by scores, if a vector (assumed to be of correct length)
map = colormap; blnFull = 0;

if (length(varargin)>1)
    scores=varargin{2};
    [r c] = size(scores);
    if c > 1
        % then full colormap passed...
        map = scores;
        blnFull = 1;
    else
        if min(scores)==0
            map=[.5 0 .5; map];
        end;
        colors = size(map,1);
        R=scores-min(scores)+1e-10; %Create R vector of scores > 0
        Rcolor=colors*R/max(R); %normalize scores/colors to number map elements
        Ucolor=unique(Rcolor);
        
        for j=1:length(Ucolor)
            idx(j)=find(Rcolor==Ucolor(j),1);
        end
        % keyboard
    end
    
end

if (length(varargin)>2) & ~isempty(varargin{3})
    shapes=varargin{3};
else
    shape='.';
    shapes=repmat(shape,length(W),1); 
end

if (length(varargin)>3)
    msize = varargin{4};
else
    msize = 10;
end
%
x=xy(:,1); y=xy(:,2);
edges=find(W);

We=[W(edges),edges];
sortWe=sortrows(We);


% This is for a Weighted Network or an unweighted network
str=(W/max(max(W))).^alpha;

% This is for Making the edges random strengths
% str=rand(size(W));


hold on
% if (length(varargin)>1)
% 
%     for i=idx,
%        
%         h=plot(x(i),y(i),shapes(i,:),'markersize',msize);
%         set(h,'Color',map(ceil(Rcolor(i)),:));            
%     end
% end
N=length(W);
for ie=sortWe(:,2)',
    i=mod(ie-1,N)+1;
    j=floor((ie-1)/N)+1;
    
    if (j>i)
        h=plot(x([i,j]),y([i,j])); hold on;
        % set(h,'color',str(i,j)*ones(1,3))
        set(h,'color',ones(1,3)- str(i,j))
        set(h,'LineWidth',.001);
        %set(h,'color',[.2 .2 .2]);
        
    end
end
if (length(varargin)>1)
    for i=1:N,
        if shapes(i,:)=='.'
            h=plot(x(i),y(i),shapes(i,:),'markersize',msize);
            if blnFull
                set(h,'MarkerFaceColor',map(i,:));
                set(h,'Color',map(i,:));
            else
                set(h,'MarkerFaceColor',map(ceil(Rcolor(i)),:));
                set(h,'Color',map(ceil(Rcolor(i)),:));

            end
        end
    end
    for i=1:N
        if shapes(i,:)~='.'
            h=plot(x(i),y(i),shapes(i,:),'markersize',msize);
            if blnFull
                set(h,'MarkerFaceColor',map(i,:));
                set(h,'Color',map(i,:));
            else
                set(h,'MarkerFaceColor',map(ceil(Rcolor(i)),:));
                set(h,'Color',map(ceil(Rcolor(i)),:));

            end

%             set(h,'MarkerFaceColor',[.7 .7 .7]);
%             set(h,'Color',map(ceil(Rcolor(i)),:));
        end
    end
else
    plot(x,y,'b.','markersize',msize)
end

%

axis equal
hold off
