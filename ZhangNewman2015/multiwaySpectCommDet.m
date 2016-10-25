function [bestNoOfGroups, bestPartition] = multiwaySpectCommDet(varargin)
% Implementation of: X Zhang, MEJ Newman (2015), Multiway spectral community
% detection in networks, Phys Rev E, (92)052808. This function tests for up
% to 50 possible groups and the best number is that of the knee of the
% approximated modularity. This implementation considers that, for a given
% tested number of groups, membership convergence is reached when 97 % or
% more vertices stopped swapping group membership.
% 
% Syntax...
% ... if ground truth is known:
% [NoGroups, partition] = multiwaySpectCommDet(Adj, groupMemb)
%
% ... if there's no known ground truth:
% [NoGroups, partition] = multiwaySpectCommDet(Adj)
%
% Where:
% Adj: binary (symmetric) adjacency matrix. The algorithm assumes an  
%      undirected, unweighted network.
% groupMemb: vector of indices of the groups to which each vertex belongs.
%      Its length must be the same as that of any of the dimensions of Adj.
% NoGroups: best number of groups detected.
% partition: vector of indices, up to NoGroups, reporting the best
%      partition found.
% 
% Ver. 1.0, Javier Caballero, 24-Oct-2016

%% SET UP DATA
% get ground-truth membership vector if given
A_adjMat = varargin{1};
if length(varargin) > 1
    g_iGroupsVerts = varargin{2};
end

% remove any auto-connections
if sum(diag(A_adjMat)) > 0
     warning(['Removed ' num2str(sum(diag(A_adjMat)))...
         ' self-edges from the '...
         'adjacency matrix.'])
    A_adjMat = abs(triu(A_adjMat) - tril(A_adjMat)); 
end

% check if we know the ground truth
if exist('g_iGroupsVerts', 'var') == 0
    % indicator that there's no ground truth available
    indGroundTruth = 0;
else
    % there is
    indGroundTruth = 1;
end


%% NETWORK STATS
% degrees of every vertex
d_degreesVerts = sum(A_adjMat, 1);

% unconnected vertices
unconnectVerts = d_degreesVerts == 0;

% remove unconnected vertices
if sum(unconnectVerts) > 0
    % modify adjacency matrix
    A_adjMat(unconnectVerts, :) = [];
    A_adjMat(:, unconnectVerts) = [];
    % if we know the ground truth
    if indGroundTruth == 1
        % modify ground truth group membership list
        g_iGroupsVerts(unconnectVerts)  = [];
    end
    % degrees of every vertex again
    d_degreesVerts = sum(A_adjMat, 1);
    % give heads-up
    warning([num2str(sum(unconnectVerts)) ' of the '...
        num2str(size(A_adjMat, 1) + sum(unconnectVerts))...
        ' original vertices were unconnected '...
        'and were ignored.'])
end

% number of vertices
n_NoVerts = size(A_adjMat, 1);

% total of edges in the network
m_totalEdges = sum(sum(triu(A_adjMat)));

% if we know the ground truth
if indGroundTruth == 1
    % starting up summation term for modularity
    sumTermQ = 0;
end

for countRow = 1:n_NoVerts
    for countCol = 1:n_NoVerts
        % modularity matrix
        B_modMat(countRow, countCol) = A_adjMat(countRow, countCol) -...
            ((d_degreesVerts(countRow)*d_degreesVerts(countCol))/...
            (2*m_totalEdges));

        % if we know the ground truth
        if indGroundTruth == 1
            % summation term for modularity
            sumTermQ = sumTermQ + B_modMat(countRow, countCol)*...
                kronDelta(g_iGroupsVerts(countRow),...
                g_iGroupsVerts(countCol));
        end
    end
end

% if we know the ground truth
if indGroundTruth == 1
    % modularity
    Q_mod = (1/(2*m_totalEdges))*sumTermQ;
end


%% VECTOR PARTITIONING ALGORITHM
% maxiumum number of iterations
maxIter = 50;

% startup values
% current number of communities
k_currentNoGroups = 1;
% current modularity
Q_currentMod = 0;
% if we know the ground truth
if indGroundTruth == 1
    % success rate vs ground truth
    normMutInfo = 0;
end
% maximum number of groups tested for
maxGroups = 50;
% pre-allocate memory for group indices matrix
g_iGroupsCurrent = nan(n_NoVerts, maxGroups);

% eigen-values and -vectors of the modularity matrix
[U_eigVects, lambda_eigVal] = eig(B_modMat, 'vector');

% flip them for decreasing eigenvalues and corresponding eigenvectors
lambda_eigVal = flipud(lambda_eigVal);
U_eigVects = fliplr(U_eigVects);



% while the maximum number of groups has not been reached
while k_currentNoGroups < maxGroups
    % current number of groups (min 2)
    k_currentNoGroups = k_currentNoGroups + 1;
    
    % rank of the analysis [k-1, n)
    % NOTE: larger p is slower but better approximates modularity value
    p_rankOfAnalys = k_currentNoGroups - 1;
    
    % (re)start stuff
    % iteration number
    NoIterations = 0;
    % indicator of group change, 0: no, 1: yes
    indGroupChange = ones(n_NoVerts, 1);
    % vertex vectors
    r_vertVect = [];  
    
    % vertex vectors, dim 1: element up to p_rankOfAnalys, dim 2: vertex 
    for countVert = 1:n_NoVerts
        for countP = 1:p_rankOfAnalys
            r_vertVect(countP, countVert) =...
                sqrt(lambda_eigVal(countP))*...
                U_eigVects(countVert, countP);
        end
    end
    
    % initial assumed group indices (random uniform)
    g_iGroupsCurrent(:, k_currentNoGroups) =...
        unidrnd(k_currentNoGroups, n_NoVerts, 1);
    
    %  while there are changes there's no convergence and the maximum 
    % number of iterations has not been reached, find the best fitting 
    % group for every vertex
    while sum(indGroupChange) > 0
        % one more iteration
        NoIterations = NoIterations + 1;
        % group vectors
        Rs_groupVect = {};
        inProd = [];
        
        for countGroup = 1:k_currentNoGroups
            % vertices in group countGroup
            iVertCurrentGroup =...
                find(g_iGroupsCurrent(:, k_currentNoGroups) == countGroup);
            % group vectors
            Rs_groupVect{countGroup} =...
                sum(r_vertVect(:, iVertCurrentGroup), 2);
            
            % inner/dot products, dim1: vertex, dim2: group
            for countVert = 1:n_NoVerts
                % if vertex countVert currently within the group countGroup
                if sum((iVertCurrentGroup - countVert) == 0) > 0% about 20X faster than ismember(countVert, iVertCurrentGroup)
                    inProd(countVert, countGroup) =...
                        (Rs_groupVect{countGroup} -...
                        r_vertVect(:, countVert))' *...
                        r_vertVect(:, countVert);
                else
                    inProd(countVert, countGroup) =...
                        Rs_groupVect{countGroup}'*...
                        r_vertVect(:, countVert);
                end
            end
        end
                
        % assign vertices to their highest inner product group
        for countVert = 1:n_NoVerts
            indMaxInProd =...
                find(inProd(countVert, :) == max(inProd(countVert, :)));
            % if the vertex is in the group(s) of its max inner product
            if sum(g_iGroupsCurrent(countVert, k_currentNoGroups) ==...
                    indMaxInProd) > 0
                % note that there was no change and don't switch membership
                indGroupChange(countVert) = 0;
            else
                % if there are two groups with equal maximum inner products...
                if numel(indMaxInProd) > 1
                    % ... pick one at random
                    indMaxInProd =...
                        indMaxInProd(unidrnd(numel(indMaxInProd)));
                end
                % put vertex in the group of its max inner product
                g_iGroupsCurrent(countVert, k_currentNoGroups) =...
                    indMaxInProd;
                % note that there was a change
                indGroupChange(countVert) = 1;
            end
        end
        
        % percentage of membership swapping after maximum iterations
        finalPctMembSwapping(k_currentNoGroups, 1) = ...
            round(sum(indGroupChange)/numel(indGroupChange)*1000)/10;
        
        % break if maximum iterations reached
        if maxIter == NoIterations
            % warn only if more than 3 % of nodes kept swapping their group
            % membership
            if finalPctMembSwapping(k_currentNoGroups) > 3
                warning(['When testing for ' num2str(k_currentNoGroups)...
                    ' groups, and after ' num2str(NoIterations) ' iterations, '...
                    num2str(finalPctMembSwapping(k_currentNoGroups))...
                    '% of vertices kept changing their membership. '...
                    'Breaking the loop and continuing.'])
            end
            break
        end
    end
    
    % approx modularity
    sumTermQ = 0;
    for countGroup = 1:k_currentNoGroups
        sumTermQ = sumTermQ + norm(Rs_groupVect{countGroup})^2;
    end
    Q_currentMod(k_currentNoGroups) = (1/(2*m_totalEdges))*sumTermQ;
    
    % if we know the ground truth
    if indGroundTruth == 1
        % current membership vector to test
        uniqueGroupI = unique(g_iGroupsCurrent(:, k_currentNoGroups));
        % if nodes were not assigned to all available groups...
        if numel(uniqueGroupI) < k_currentNoGroups
            % ... make a membership vector compatible with MIpartitions
            g_iGroups2Test = nan(n_NoVerts, 1);
            for count = 1:numel(uniqueGroupI)                
                g_iGroups2Test(g_iGroupsCurrent(:, k_currentNoGroups) ==...
                    uniqueGroupI(count)) = count;
            end
        else% ... use it as it is
            g_iGroups2Test = g_iGroupsCurrent(:, k_currentNoGroups);
        end
        % compute normalized mutual info
        normMutInfo(k_currentNoGroups) =...
            MIpartitions(g_iGroupsVerts, g_iGroups2Test);
    end
    
    noEffectComm(k_currentNoGroups) = numel(find(sum(cell2mat(Rs_groupVect), 1) ~= 0));
end

% best number of groups as the detected modularity knee
bestNoOfGroups = knee_pt(Q_currentMod, 1:numel(Q_currentMod));
% detected best partition, indices starting at 1, unconnected nodes'
% membership is reported as NaN
dummy = nan(n_NoVerts, 1);
uniqueGroupI = unique(g_iGroupsCurrent(:, bestNoOfGroups));
for count = 1:numel(uniqueGroupI)
    dummy(g_iGroupsCurrent(:, bestNoOfGroups) ==...
        uniqueGroupI(count)) = count;
end
bestPartition = nan(numel(unconnectVerts), 1);
bestPartition(unconnectVerts == 0) = dummy;

% report results in command window
disp(' ')
disp('RESULTS')
disp('-------')
% if we know the ground truth
if indGroundTruth == 1
    disp('GROUND TRUTH')
    disp(['Modularity: ' num2str(Q_mod)])
    disp(['No. of groups: ' num2str(max(g_iGroupsVerts))])
    disp(' ')
end
disp('FOUND')
disp(['Max approx. modularity reached: '...
    num2str(max(Q_currentMod))])
disp(['Best No. of groups: ' num2str(bestNoOfGroups)])
disp(['Approx. modularity at best No of groups: '...
    num2str(Q_currentMod(bestNoOfGroups))])
% if we know the ground truth
if indGroundTruth == 1
    disp(['Norm. mutual info. at best No. of groups: '...
        num2str(normMutInfo(bestNoOfGroups))])
end























