% unit-test consensus clustering of network functions
clear all; close all
addpath ../

%% make data
% uniformly random weighted network
n = 100;
W = rand(100); % uniformly random weight matrix
W = W+W'/2;    % undirected 
W(eye(n)==1) = 0; % remove self-connections
Data{1} = W;

% block diagoonal
Wblock = ones(n/2); 
W = zeros(n);
W(1:n/2,1:n/2) = Wblock;
W(n/2+1:end,n/2+1:end) = Wblock;
W(eye(n)==1) = 0; % remove self-connections
Data{2} = W;

%% parameter options
dims = {'all','scale'};
Treps = [0,1,10];

%% test data 
for iD = 1:numel(Data)
    % get modularity: do this with configuration model to get appropriate data
    B = Data{iD} - expectedA(Data{iD});  % also catch expectedA errors here
    
    % get eigenvalues & eigenvectors 
    [V,egs] = eig(B,'vector');
    [egs,ix] = sort(egs,'descend');    % sort into descending order
    V = V(:,ix);                       % ditto the eigenvectors 
 
    % get embedding dimensions
    ixRetain = 1:3;
    % ixRetain = find(egs > 0);        % gives huge numbers for noise
    EmbedD = V(:,ixRetain);
    M = numel(ixRetain)+1;  % how many groups
    
    %% test k-means sweep
    Mlist = [0,1,M,M+5];  % test passing of no groups, 1 group, matching groups, and more than dimensions
    ErrCtr = 1;
    for iM = 1:numel(Mlist)
        for iT = 1:numel(Treps)
            for iOpt = 1:numel(dims)
                try
                    allgrps = kmeansSweep(EmbedD,Mlist(iM),Treps(iT),dims{iOpt});
                catch ME
                    msg = sprintf(['\n Data' num2str(iD) ', K-means Error' num2str(ErrCtr) ': \n' ME.message ' \n Thrown by: M=' num2str(Mlist(iM)) ', Treps = ' num2str(Treps(iT)) ', dims = ' dims{iOpt}]); 
                    disp(msg)
                    allgrps = [];
                    ErrCtr = ErrCtr+1;
                end
                % sanity checks
                if ~isempty(allgrps)
                    [r,c] = size(allgrps);
                    if c ~= Mlist(iM)*Treps(iT) || r ~= n;
                        error('K-means sweep output is the wrong size')
                    end
                end
                if Mlist(iM) == M && Treps(iM) == max(Treps) && any(strfind(dims{iOpt},'all'))
                    Clustering = allgrps;  % use this for all next tests
                end
            end
        end
    end
    
    %% Q computation
    m = sum(sum(Data{iD}))/2;  % sum of unique weights
    try
        Q = computeQ(Clustering(:,1),B,m);
    catch ME
        msg = sprintf(['\n Data' num2str(iD) ', Q compute Error: \n' ME.message ]); 
        disp(msg)
    end
    
    %% 
end
