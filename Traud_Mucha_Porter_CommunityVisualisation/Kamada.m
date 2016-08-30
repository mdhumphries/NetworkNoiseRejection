function xynew = Kamada(A,epsilon)
% This code follows the algorithm in the article written by Tomihisa Kamada
% and Satoru Kawai entitled "An Algorithm For Drawing General Undirected
% Graphs."
%
% Inputs:
%
% A is the adjacency matrix of your network
%
% epsilon is the accuracy at which the algorithm acts, .01 works well.
%
% Outputs:
%
% xynew is the matrix of xy coordinates for each node.
%
% Calling this Function:
%
% xynew = Kamada(A,epsilon)
%
% Last Modified by ALT May 14, 2009
maxiter=5;
A = sparse(A);
A(find(A))=1./A(find(A));
%Made to calculate Kamada Kawai values from AN Algorithm for Drawing
%General Undirected Graphs
n=length(A);





% compute d_ij
d_ij=floyd_warshall_all_sp(A);

% compute l_ij
%l=l_ij(L0,d_ij);

l=d_ij;
% compute k_ij
k=k_ij(d_ij);
% initialize p's

xy=rand(n,2)-.5;
dE_dx=zeros(n,1);
dE_dy=zeros(n,1);
for m=1:n
    [dE_dx(m) dE_dy(m)] = Eq78(xy,k,l,m);
end
del=(((dE_dx.^2)+(dE_dy.^2)).^(.5));

xynew=xy;
counterm=0;

mtrac = [];
% While loop of max_i(del_i)>Epsilon
while(max(del)>epsilon)

    % let p_m be the particle satisfying del_m = max_i(del_i)
    m=find(max(del)==del);
    if length(m)~=1
        m = m(1);
    end
    mtrac(end + 1) = m;
    l2 = length(mtrac);
   if l2>=6
       vector = m*ones(1,6);
       l6 = mtrac(l2-5:end);
       counterm = sum(vector==l6);
       

    if counterm>=6
        xym = xynew(m,:);
        xynew(m,:) = (xynew(m,:)-rand(1,2)-.5);
        [dE_dx(m) dE_dy(m)]=Eq78(xynew, k,l,m);


        index=1:n;
        index=find(index~=m);
        for j=index
            [dE_dx(j) dE_dy(j)]=recalcEq78(dE_dx(j), dE_dy(j),xym,xynew, j,m,l,k);
        end
        del=(((dE_dx.^2)+(dE_dy.^2)).^(.5));

        m=find(max(del)==del);
        if length(m)~=1
            m = m(1);
        end

    end 
    end
    % While loop of del_m > Epsilon
    xy=xynew; %
    niter=0;
    while((del(m)>epsilon)&(niter<maxiter))
        niter=niter+1;
    %for iiter=1:2,
        %Calculate dE_dx and dE-dy for this particular m
        
        %Calculate the second derivatives for this particular m
        [d2E_dx2 d2E_dxdy d2E_dydx d2E_dy2] = Eq1316(m, k, l, xynew);
        %Use the derivatives and second derivatives to solve for del_x and del_y
        
        x=[d2E_dx2, d2E_dxdy; d2E_dydx, d2E_dy2]\[-dE_dx(m);-dE_dy(m)];
        del_x=x(1); del_y=x(2);
        % x_m = x_m + del_x
        xynew(m,1)=xynew(m,1)+del_x;
        % y_m = y_m + del_y
        xynew(m,2)=xynew(m,2)+del_y;
        %Should be only recalculating one delta
        
            %[dE_dx(m) dE_dy(m)] = Eq78(xynew,k,l,m);
        [dE_dx(m) dE_dy(m)]=Eq78(xynew, k,l,m);
        del(m)=(((dE_dx(m)^2)+(dE_dy(m)^2))^(.5));
        
    end
    xym=xy(m,:);
    %recalcEq78
   index=1:n;
   index=find(index~=m);
    for j=index
        [dE_dx(j) dE_dy(j)]=recalcEq78(dE_dx(j), dE_dy(j),xym,xynew, j,m,l,k);
    end
    [dE_dx(m) dE_dy(m)]=Eq78(xynew,k,l,m);
    del=(((dE_dx.^2)+(dE_dy.^2)).^(.5));
    
    
end


end


%--------------------------------------------------------------------------
function [k] = k_ij(d_ij)
%Calculates k_ij for the Kamada Kawai code from d_ij

% Set K equal to 1; actual value does not effect the outcome
K = 1;
k = K./d_ij.^2;
% make the diagonal of k equal to zero (there isn't a spring bewteen a node
% and itself
for i = 1 : length(k)
    k(i,i) = 0;
end

end
%--------------------------------------------------------------------------
function [dE_dx dE_dy] = Eq78(xy,k,l,m,i,xynew)

% This function computes formula (7) and (8) from the Kamada Kawai paper
%       for loop on the outside
n=length(xy);
if(nargin==6)
    index=i;
    xi=xynew(index,1);
    yi=xynew(index,2);
else
    index = 1:n;
    index = index(find(index~=m));
    xi = xy(index,1);
    yi = xy(index,2);
end

xm=xy(m,1);
ym=xy(m,2);


lmi = l(m,index)';
kmi = k(m,index)';

a1 = (xm - xi);
a2 = (ym - yi);
frac=lmi./(((a1.^2)+(a2.^2)).^.5);
dE_dx = sum(kmi.*(a1-(a1.*frac)));
dE_dy = sum(kmi.*(a2-(a2.*frac)));
end
%--------------------------------------------------------------------------
function [dE_dx dE_dy]=recalcEq78(dE_dx, dE_dy, xym, xynew, m, oldm,l,k)


% Subtract old contribution
dE_dxold=(k(m,oldm)*((xynew(m,1)-xym(1))-((l(m,oldm)*(xynew(m,1)-xym(1)))/sqrt(((xynew(m,1)-xym(1))^2)+((xynew(m,2)-xym(2))^2)))));
dE_dyold=(k(m,oldm)*((xynew(m,2)-xym(2))-((l(m,oldm)*(xynew(m,2)-xym(2)))/sqrt(((xynew(m,1)-xym(1))^2)+((xynew(m,2)-xym(2))^2)))));
dE_dx=dE_dx-dE_dxold;
dE_dy=dE_dy-dE_dyold;
% Add in new contribution
xym=xynew(oldm,:);
dE_dxnew=(k(m,oldm)*((xynew(m,1)-xym(1))-((l(m,oldm)*(xynew(m,1)-xym(1)))/sqrt(((xynew(m,1)-xym(1))^2)+((xynew(m,2)-xym(2))^2)))));
dE_dynew=(k(m,oldm)*((xynew(m,2)-xym(2))-((l(m,oldm)*(xynew(m,2)-xym(2)))/sqrt(((xynew(m,1)-xym(1))^2)+((xynew(m,2)-xym(2))^2)))));
dE_dx=dE_dx+dE_dxnew;
dE_dy=dE_dy+dE_dynew;
end
%--------------------------------------------------------------------------

function [d2E_dx2 d2E_dxdy d2E_dydx d2E_dy2] = Eq1316(m, k, l, xy)

% This function calculates the second derivatives needed to find delta_x
% and delta_y in the kamada kawai algorithm

n = length(xy);

index = 1:n;
index = index(find(index~=m));

xi = xy(index,1);
yi = xy(index,2);
xm=xy(m,1);
ym=xy(m,2);
lmi = l(m,index)';
kmi = k(m,index)';

a1 = (xm - xi);
a2 = (ym - yi);
frac=lmi./(((a1.^2)+(a2.^2)).^1.5);
d2E_dx2 = sum(kmi.*(1-(frac.*(a2.^2))));
d2E_dxdy = sum(kmi.*(((a1).*(a2)).*frac));
d2E_dydx = d2E_dxdy;
d2E_dy2 = sum(kmi.*(1-(frac.*(a1.^2))));
end
%--------------------------------------------------------------------------

