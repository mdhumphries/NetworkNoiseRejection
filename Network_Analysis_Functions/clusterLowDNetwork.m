function Results = clusterLowDNetwork(W,P,N,M,nreps,group_membership)

% CLUSTERLOWDNETWORK clusters low-dimensional projection of network
% R = CLUSTERLOWDNETWORK(W,P,N,M,nreps,G) clusters the (n*n) network W using
% its projection into a low-dimensonal space, between N and M dimensions
% (often N==M). Also needs:
%       P : expected null model used to define the low-D space
%       nreps : number of repeats of k-means at each number of dimensions
%       G : ground-truth clustering as an n-length array
% 
% 9/7/2018 : initial version
% Mark Humphries 
 

 if N > 1     % if dimensions exist to cluster in...
    [QmaxCluster,Qmax,ConsCluster,ConsQ,~] = ...
              ConsensusCommunityDetect(W,P,N,M,nreps);
    % quality of estimation of retained communities
    if ~isempty(QmaxCluster)
        [~,Results(iB).normVIQmaxSpectra(iP)]=VIpartitions(QmaxCluster,group_membership);
        Results(iB).nGrpsQmaxSpectra(iP) = max(QmaxCluster);
    else
        Results(iB).normVIQmaxSpectra(iP)= 0;
        Results(iB).nGrpsQmaxSpectra(iP) = nan;                
    end

    if ~isempty(ConsCluster)
        [~,Results(iB).normVIConsensusSpectra(iP)]=VIpartitions(ConsCluster,group_membership);
        Results(iB).nGrpsConsensusSpectra(iP) = max(ConsCluster);
    else
        Results(iB).normVIConsensusSpectra(iP)= 0;
        Results(iB).nGrpsConsensusSpectra(iP) = nan;                
    end
else
    Results(iB).normVIQmaxSpectra(iP)=0;
    Results(iB).normVIConsensusSpectra(iP)=0;
    Results(iB).nGrpsQmaxSpectra(iP) = 0;
    Results(iB).nGrpsConsensusSpectra(iP) = 0;
end
