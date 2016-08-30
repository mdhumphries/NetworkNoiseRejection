function [xynew, XY, xyc] = fruc_reinC(A, gn, epsilon, seed, factor1, factor2)

% This function uses the Fruchterman-Reingold Algorithm to find the optimal
% node placement for a given network with respect to Communities.
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
%[xynew, XY, xyc] = fruc_reinC(A, gn, epsilon, seed, factor1, factor2)
%
% Last Modified by ALT May 14, 2009
%
%Modified by Lucas on November 14, 2011

%deal with row/column vector
if (size(gn,2)==1)
    gn=gn';
end


CAM=commAdjMatrixSparse(gn,A);
gnu=unique(gn);
num_com=length(unique(gn));
rand('twister',seed);

nodes_percom = zeros(num_com,1);
for i = 1:length(gnu)
    nodes_percom(i) = length(find(gnu(i)==gn));
end
P = (CAM./(nodes_percom*nodes_percom'));
P=P./max(max(P));
P=10*(P-diag(diag(P)));


%Calculate Kamada for Community Placement

XY=fruc_rein(P,epsilon);
XY=factor1*XY;


%Calculate minA for all communities

xyc=cell(1,num_com);

% Find the distances between communities
tic
for j=1:num_com
    minA=MinAdjMatrix(A,gnu(j), gn);
    
    if nodes_percom(j)==1
        xyc{j}=[0 0];
    else
    xyc{j}=fruc_rein(minA,epsilon);
    end
    xyca=xyc{j};
    %Make xyc centered around zero
    xyca(:,1)=xyca(:,1)-mean(xyca(:,1));
    xyca(:,2)=xyca(:,2)-mean(xyca(:,2));
    if max(max(abs(xyca)))>2
        xyca=factor2*(xyca/max(max(abs(xyca))));
    end
    xyc{j}=5*xyca;

end
t2=toc
tic
% Make KamadaR to calculate the rotation or flip of each community
XY=XY*.01*length(gn);

xynew=frun_reinR(A,xyc,XY,gn,epsilon);

t3=toc



end
%--------------------------------------------------------------------------
function [xy] = fruc_rein(A,eps)
% This function uses the Fruchterman-Reingold Algorithm to find the optimal
% node placement for a given network.
% A is the adjacency matrix for the network


n = length(A);
W = 3;
t=50;
L = 3;
iter=0;
maxiter=length(A);
% Initialize the point's
xy=rand(n,2)-.5;
k = .75*sqrt((W*L)/n);
displacement = zeros(n,2);
displacement(1,1)=1;
z = 0;
while (((max(max(abs(displacement))))>=eps)&&(iter<maxiter))
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
        displacement(row(l),:) = displacement(row(l),:)- A(row(l),col(l))*(((d/norm(d))*fa(norm(d),k))*.5);
        displacement(col(l),:) = displacement(col(l),:)+ A(col(l),row(l))*(((d/norm(d))*fa(norm(d),k))*.5);
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
    iter=iter+1;
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
            %comm1=comm1(find(comm1));
            %comm2=comm2(find(comm2));
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
    %Is there a more efficient way to specify this in MATLAB?
    communities=groups;
    for ig=1:length(groupnumbers),
        indx=find(groups==groupnumbers(ig));
        communities(indx)=commnumbers(ig);
    end
end
end

%--------------------------------------------------------------------------

function [minA] = MinAdjMatrix(A, gn, gnlist)
% This function will make a adjacentcy matrix for a particular community
% that is in the original adjacentcy matrix (A) and indicated by group number
% (gn) from a group list (gnlist)

index = find(gn==gnlist);
minA = A(index,index);

end
%--------------------------------------------------------------------------
function xy=frun_reinR(A,xyc,XY,gn,epsilon)


n = length(A);
W = 3;

L = 3;
% Initialize the point's

k = .75*sqrt((W*L)/n);


% Make xy coordinates
gnu=unique(gn);
xy=zeros(length(gn),2);
idx=cell(1,length(gnu));
n=length(gn);
displacement = zeros(n,2);
dispc(1)=1+epsilon;
maxiter=0;
while((max(dispc)>epsilon)&& maxiter<(length(gnu)))
    maxiter=maxiter+1;
    for i=1:length(gnu)
        idx{i}=find(gnu(i)==gn);

        xy(idx{i},1)=XY(i,1)+xyc{i}(:,1);
        xy(idx{i},2)=XY(i,2)+xyc{i}(:,2);
    end;
    %Calculate Fa
    [row,col] = find(A);
    displacement = zeros(n,2);
    for l = 1:length(row)
        d = xy(row(l),:)-xy(col(l),:);
        displacement(row(l),:) = displacement(row(l),:)-(((d/norm(d))*fa(norm(d),k))*.5);
        displacement(col(l),:) = displacement(col(l),:)+(((d/norm(d))*fa(norm(d),k))*.5);
    end
    
    
    
    % Calculate Fr
    for i = 1:n
        indx=1:n;
        indx=find(indx~=i);
        xyi=repmat(xy(i,:),n-1,1);
        d=xyi-xy(indx,:);
        normd=((d(:,1).^2)+(d(:,2).^2)).^.5;
        displacement(i,:)=displacement(i,:)+sum([(d(:,1)./normd).*fr(normd,k) (d(:,2)./normd).*fr(normd,k)]);
    end
   
    % displacement is now the difference between Fa and Fr
    dispc=zeros(length(gnu),1);
    disple=((displacement(:,1).^2)+(displacement(:,2).^2)).^.5;
    for i=1:length(gnu)
        dispc(i)=sum(disple(idx{i}));
    end
    
    
    i=find(max(dispc)==dispc);
    comm=i;
    try
    comms(maxiter)=comm;
    catch
        keyboard
    end
    if maxiter>3
        commse=repmat(comms(end),1,3);

        if (sum(commse==comms(end-2:end))==3)
            i=ceil(rand()*length(dispc));

            comm=i;
            comms(maxiter)=comm;

        end
    end
  
    
    theta=0;

    [thetaR, FDiffR] = fsolve(@(theta)thetafunc(theta,A,xyc,gn,XY,idx,comm,0),theta);

    [thetaF, FDiffF] = fsolve(@(theta)thetafunc(theta,A,xyc,gn,XY,idx,comm,1),theta);
   
    Cdiff=[FDiffR,FDiffF];
    Ctheta=[thetaR,thetaF];

    newtheta=Ctheta(find(min(Cdiff)==Cdiff));
    if length(newtheta==2)
        newtheta=newtheta(1);
    end
    
    if ((find(min(Cdiff)==Cdiff))==2)
        xyc{comm}=xyc{comm}*[-cos(newtheta), -sin(newtheta); -sin(newtheta), cos(newtheta)];
        xy(idx{comm},1)=XY(comm,1)+xyc{comm}(:,1);
        xy(idx{comm},2)=XY(comm,2)+xyc{comm}(:,2);
        dispc(comm)=FDiffF;
    else
        try
        xyc{comm}=xyc{comm}*[cos(newtheta), sin(newtheta); -sin(newtheta), cos(newtheta)];
        catch
            keyboard
        end
        xy(idx{comm},1)=XY(comm,1)+xyc{comm}(:,1);
        xy(idx{comm},2)=XY(comm,2)+xyc{comm}(:,2);
        dispc(comm)=FDiffR;
    end

end


end

%--------------------------------------------------------------------------

function [Fdiff]=thetafunc(theta,A,xyc,gn,XY,idx,com,flip)

if flip
    xyc{com}(:,1)=-(xyc{com}(:,1).*cos(theta))-(xyc{com}(:,2).*sin(theta));
    xyc{com}(:,2) = -(xyc{com}(:,1).*sin(theta))+(xyc{com}(:,2).*cos(theta));

else
    xyc{com}(:,1)=(xyc{com}(:,1).*cos(theta))-(xyc{com}(:,2).*sin(theta));
    xyc{com}(:,2) = (xyc{com}(:,1).*sin(theta))+(xyc{com}(:,2).*cos(theta));

end
%calculate Fdiff
gnu=unique(gn);
for i=1:length(gnu)
    xy(idx{i},1)=XY(i,1)+xyc{i}(:,1);
    xy(idx{i},2)=XY(i,2)+xyc{i}(:,2);
end
W = 3;

L = 3;
% Initialize the point's
n=length(A);
k = .75*sqrt((W*L)/n);
%Calculate Fa
[row,col] = find(A);
id=idx{com};




displacement = zeros(n,2);


for l = id
   


    row2=row(find(row==l));
    col2=col(find(row==l));
   
    for i=1:length(row2)
    d = xy(row2(i),:)-xy(col2(i),:);
    displacement(row2(i),:) = displacement(row2(i),:)-(((d/norm(d))*fa(norm(d),k))*.5);
    displacement(col2(i),:) = displacement(col2(i),:)+(((d/norm(d))*fa(norm(d),k))*.5);
    end
   
    col3=col(find(col==l));
    row3=row(find(col==l));
    for i=1:length(row2)
    d = xy(row3(i),:)-xy(col3(i),:);
    displacement(row3(i),:) = displacement(row3(i),:)-(((d/norm(d))*fa(norm(d),k))*.5);
    displacement(col3(i),:) = displacement(col3(i),:)+(((d/norm(d))*fa(norm(d),k))*.5);
    end
   
end

% Calculate Fr


for i = (idx{com})
    indx=1:n;
    indx=find(indx~=i);
    xyi=repmat(xy(i,:),n-1,1);
    d=xyi-xy(indx,:);
    normd=((d(:,1).^2)+(d(:,2).^2)).^.5;
    displacement(i,:)=displacement(i,:)+sum([(d(:,1)./normd).*fr(normd,k) (d(:,2)./normd).*fr(normd,k)]);
end

disple=((displacement(:,1).^2)+(displacement(:,2).^2)).^.5;

dispc=sum(disple(idx{com}));

Fdiff=dispc;
end

