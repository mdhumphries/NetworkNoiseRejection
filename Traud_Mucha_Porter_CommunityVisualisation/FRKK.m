function [xynew, XY, xyc]=FRKK(A,gn,epsilon, seed, factor1, factor2)

% This function uses the Fruchterman-Reingold Algorithm to find the optimal
% node placement for a given network with respect to Communities.  But then
% uses the Kamada Kawai algorithm to place the nodes within each community.
%
% Inputs:
%
% A is the adjacency matrix for the network
%
% gn are a vector of group numbers, i.e. every node that has the number 1
% is in group 1.
%
% episilon is the accuracy at which you would like the algorithm to act,
% .01 tends to work well.
%
% seed is the number used in the random number generator for placing the
% communities, any four digit number will work.
%
%factor1 is the number to multiply the community coordinates by, 2 tends to
%work well, but you may want the communities farther apart.
%
%factor2 is the number to multiply the node coordinates by within the
%community, if there are more than five nodes within a particular community
%it is helpful to normalize the coordinates and then multiply them by a
%factor to redistribute them.  4 tends to work in most cases but you may
%want to change it based on your data.
%
%Outputs:
%
%xynew are the vectors of the xy coordinates for each node in your network, this
%is one of the inputs for graphplot2d.m
%
% XY are the vectors of xy coordinates for the communties, this is one of
% the inputs for drawForceCPie.m
%
%xyc is the cellarray of xy coordinates for the nodes within each community
%that are centered at zero.
%
%Calling this Function:
%
%[xynew, XY, xyc]=FRKK(A,gn,epsilon)
%               or
%[xynew, XY, xyc]=FRKK(A,gn,epsilon, seed, factor1, factor2)
%
%
% Last Modified by ALT May 14, 2009

% Distances are written as one over the weight for a weighted matrix (the
% more weight representing a link, the closer in distance they should be)
%
%Last Modified by Lucas November 15, 2011 to fix bug on line 362 

%check wether group vector is row or column
if (size(gn,2)==1)
    gn=gn';
end

if nargin==3
    seed=1327;
    factor1=2;
    factor2=4;
end
CAM=commAdjMatrixSparse(gn,A); % the total number of edges between each pair of communities

gnu=unique(gn);
num_com=length(unique(gn));
rand('twister',seed);

nodes_percom = zeros(num_com,1);
for i = 1:length(gnu)
    nodes_percom(i) = length(find(gnu(i)==gn));
end
P =(CAM./((nodes_percom*nodes_percom')));
%P=P./(max(max(P))+1);
P=(P-diag(diag(P))); % the proportion of all possible edges between each pair of communities, 
                     % minus self-connections (which should be 0 anyway at this point)

%Calculate Kamada for Community Placement

XY=fruc_rein(P,epsilon);
XY=factor1*XY;
m=max(nodes_percom);
%Calculate minA for all communities

xyc=cell(1,num_com);


% Find the distances between communities

for j=1:num_com
    minA=MinAdjMatrix(A,gnu(j), gn);
    xyc{j}=Kamada(minA,epsilon);
    xyca=xyc{j};
    %Make xyc centered around zero
    xyca(:,1)=xyca(:,1)-mean(xyca(:,1));
    xyca(:,2)=xyca(:,2)-mean(xyca(:,2));
    alpha=(length(xyca)/m)*factor2;
    if length(xyca)>5
        xyc{j}=alpha*(xyca/max(max(abs(xyca))));
    end
end


% Make KamadaR to calculatethe rotation or flip of each community

xynew=KamadaR(A,xyc,XY,gn,epsilon);




end
%--------------------------------------------------------------------------
function [xy] = fruc_rein(A,eps)
% This function uses the Fruchterman-Reingold Algorithm to find the optimal
% node placement for a given network.
% A is the adjacency matrix for the network



n = length(A);
W = 5;
t=50;
L = 5;
% Initialize the point's
xy=rand(n,2)-.5;
k = .75*sqrt((W*L)/n);
displacement = zeros(n,2);
displacement(1,1)=1;
z = 0;
while ((max(max(abs(displacement))))>=eps)
    [row,col] = find(A);
    z=z+1;
    
    m(z)=max(max(abs(displacement)));
    l2=length(m);
    if l2>=6
       v=diff(m);
       

    if abs(v(end))<=.00000001
        break
    end 
    end
    displacement = zeros(n,2);
    
    for l = 1:length(row)
        d = xy(row(l),:)-xy(col(l),:);
        displacement(row(l),:) = displacement(row(l),:)-A(row(l),col(l))*(((d/norm(d))*fa(norm(d),k))*.5);
        displacement(col(l),:) = displacement(col(l),:)+A(col(l),row(l))*(((d/norm(d))*fa(norm(d),k))*.5);
    end
    for i = 1:n
        idx=1:n;
        idx=find(idx~=i);
        xyi=repmat(xy(i,:),n-1,1);
        d=xyi-xy(idx,:);
        normd=((d(:,1).^2)+(d(:,2).^2)).^.5;
        displacement(i,:)=displacement(i,:)+sum([(d(:,1)./normd).*fr(normd,k) (d(:,2)./normd).*fr(normd,k)]);
        xy(i,:) = xy(i,:) + (displacement(i,:)/norm(displacement(i,:)))*min(norm(displacement(i,:)),t);
    end
    t=.95*t;

end

end

%--------------------------------------------------------------------------
function [f_a] = fa(distance,k)

f_a = ((distance^2)./k);
end

%--------------------------------------------------------------------------
function [f_r] = fr(distance,k)


f_r = ((k^2)./distance);
end

%--------------------------------------------------------------------------

function mat = commAdjMatrixSparse(groups, A)
% Creates a community adjacency matrix using
% groups from the output for reccurrcommsNew2Sparse, A is the adjacency matrix
% 0's on the diagonal, other elements consist of the total number of
% connections between the two communities
%
% Last Modified by ALT 20 June 2007


h=sort(groups);
g=unique(h);
d=diff(g);
f=sort(d);
z=unique(f);
cuts=size(z,2);

[communities cut]=findcommunitiesatcut(groups,cuts);
rows = max(communities);
mat=spalloc(rows,rows,2*rows);

for i = 1:rows
    for j = 1:rows
        if(i ~= j)
            comm1 = find(communities==i);
            comm2 = find(communities==j);
            mat(j, i) = sum(sum(A(comm1, comm2)));
        end
    end
end
end


%-------------------------------------------------------------------------
function [communities cut] = findcommunitiesatcut(groups,cut)
%[communities cut]=findcommunitiesatcut(groups,cut)
%
% Gives the community numbers at a requested cut or level in the groups vector,
% if the cut number is not valid the program changes it to a valid one.
% Uses a groups vector and a scalar cut number, gives communities and the cut number,
% which is needed when cut is changed.
%
%
%Last modified by ALT, 20 June 2007

%Error checking
n=unique(groups);
f=diff(n);
z=unique(f);
cutmax=length(z);
if(cut>cutmax)
    disp(['That is too many cuts! I have changed the cut number.']);
    cut=cutmax;
elseif(cut<0)
    cut=cutmax;
    disp(['Negative numbers dont work, I have changed the cut number to the max'])
end

%Identify distinct group values and number of cut levels in dendrogram:
groupnumbers=unique(groups);
differences=diff(groupnumbers);
diffnumbers=unique(differences);
cuts=length(diffnumbers);

if cut==0,
    communities=ones(size(groups));
else
    cutdiff=diffnumbers(cuts+1-cut); %NOTE THERE IS NO ERROR CHECKING HERE, ASSUMED VALID CUT NUMBER
    
   
    commnumbers=cumsum([1,diff(groupnumbers)>=cutdiff]);

    %Define communities by replacing the groupnumbers values in groups with the
    %corresponding commnumbers values, component by component.
    communities=groups;
    for ig=1:length(groupnumbers),
        indx=find(groups==groupnumbers(ig));
        communities(indx)=commnumbers(ig);
    end
end
end
%----------------------------------------------
function [xynew] =KamadaR(A,xyc,XY,gn,epsilon)



%Make full xy coordinates
gnu=unique(gn);
xy=zeros(length(gn),2);
idx=cell(1,length(gnu));
for i=1:length(gnu)
    idx{i}=find(gnu(i)==gn);
    xy(idx{i},1)=XY(i,1)+xyc{i}(:,1);
    xy(idx{i},2)=XY(i,2)+xyc{i}(:,2);
end;


d_ij=floyd_warshall_all_sp(A);
l=d_ij;
k=k_ij(d_ij);

%Calculate dE_dx and dE_dy
n=length(gn);
dE_dx=zeros(n,1);
dE_dy=zeros(n,1);
for m=1:n
    [dE_dx(m) dE_dy(m)] = Eq78(xy,k,l,m);
end
theta=0;
dE_dtheta=zeros(1,length(gnu));


for j=1:length(gnu)
    dE_dtheta(j)=thetafunc(theta,xyc{j},dE_dx,dE_dy,idx{j},0);
end

maxiter=0;
while((max(dE_dtheta)>epsilon)&& maxiter<(length(gnu)))

    %Find Community to rotate

    comm=find(max(dE_dtheta)==dE_dtheta,1);







    %Find angle to rotate it or if to flip and rotate
    [thetaR, dE_dthetaR] = fsolve(@(theta)thetafunc(theta,xyc{comm},dE_dx,dE_dy,idx{comm},0),theta);

    [thetaF, dE_dthetaF] = fsolve(@(theta)thetafunc(theta,xyc{comm},dE_dx,dE_dy,idx{comm},1),theta);
    CdE_dtheta=[dE_dthetaR,dE_dthetaF];
    Ctheta=[thetaR,thetaF];

    newtheta=Ctheta(find(min(CdE_dtheta)==CdE_dtheta));
    if ((find(min(CdE_dtheta)==CdE_dtheta))==2)
        rxy=xyc{comm}*[-cos(newtheta), -sin(newtheta); -sin(newtheta), cos(newtheta)];
        xy(idx{comm},1)=XY(comm,1)+rxy(:,1);
        xy(idx{comm},2)=XY(comm,2)+rxy(:,2);
    else
        rxy=xyc{comm}*[cos(newtheta), sin(newtheta); -sin(newtheta), cos(newtheta)];
        xy(idx{comm},1)=XY(comm,1)+rxy(:,1);
        xy(idx{comm},2)=XY(comm,2)+rxy(:,2);
    end
    for m=1:n
        [dE_dx(m) dE_dy(m)] = Eq78(xy,k,l,m);
    end
    for j=1:length(gnu)
        dE_dtheta(j)=thetafunc(theta,xyc{j},dE_dx,dE_dy,idx{j},0);
    end
    maxiter=maxiter+1;
end
xynew=xy;


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


%n=length(xy);

%changed to avoid error when passing a single point - Lucas
n=size(xy,1);
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
function [dE_dtheta]=thetafunc(theta,xyc,dE_dx,dE_dy,idx,flip)

if flip
    dx_dtheta=(xyc(:,1).*sin(theta))-(xyc(:,2).*cos(theta));
    dy_dtheta = -(xyc(:,1).*cos(theta))-(xyc(:,2).*sin(theta));

else
    dx_dtheta=-(xyc(:,1).*sin(theta))-(xyc(:,2).*cos(theta));
    dy_dtheta = (xyc(:,1).*cos(theta))-(xyc(:,2).*sin(theta));

end
dE_dtheta=dot(dE_dx(idx),dx_dtheta)+dot(dE_dy(idx),dy_dtheta);
end

%--------------------------------------------------------------------------

function xynew = Kamada(A,epsilon)
rand('twister',1313);
maxiter=5;
maxiter2=10; 

A = sparse(A);
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
counterm = 0;

mtrac = [];
% While loop of max_i(del_i)>Epsilon
iter2=0;
while((max(del)>epsilon)&&(iter2<maxiter2))

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
    while((del(m)>epsilon)&&(niter<maxiter))
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
    iter2=iter2+1;

end


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
function MinA=MinAdjMatrix(A,gnu, gn)
idx=find(gn==gnu);
MinA=A(idx,idx);




end