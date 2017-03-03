function allgrps = kmeansSweep(D,K,Treps,dims)

% KMEANSSWEEP sweep of kmeans clustering across multiple K
% C = KMEANSSWEEP(D,K,R,DIMS) performs a k-means sweep of a data matrix, given:  
% D: set of data points (nxd): n objects in d dimensions  (e.g. low-D
% projection)
% K: maximum number of groups to test (minimum assumed to be 2)
% R: repeats of each clustering at each K
% DIMS = 'all' or 'scale': 'all' uses all of D, irrespective of K; 'scale'
% uses only the i-1 set of dimensions
%
% Returns C, a m x [(K-1)xR] matrix of every k-means clustering, in order of low to high K 
%   (columns 1:K-1 are for K=2; columns K to 2*K-1 are for K=3 etc) 
%
% Mark Humphries 2/3/2017

if K < 2
    error('Specify at least K=2 groups')
end

[n,d] = size(D);
if any(strfind(dims,'scale')) && d < K-1
    error('Not enough embedding dimensions to scale to K')
end

ixNow = 1;
allgrps = zeros(n,(K-1)*Treps);
for ngrps = 2:K
    switch dims
        case 'all'
            thisVector = D;
        case 'scale'
            thisVector = D(:,ngrps-1); % subset of dimensions
        otherwise
            error('Unknown options for DIMS')
    end
    
    cpos = kmeansplus(thisVector, ngrps); % initialise centers
    for rep = 1:Treps
        % keyboard
        try            
            allgrps(:,ixNow) = kmeans(thisVector,ngrps,'Distance','sqeuclidean','Start',cpos);

        catch
            % if kmeans throws a wobbly, set to "no groups"...
            warning('kmeans wobbly')
        end
        ixNow = ixNow + 1;
    end % end k-means repeats loop
end % end groups loop