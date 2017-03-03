% unit-test consensus clustering of network functions
clear all; close all
addpath ../

egmin = 1e-2; % counts as more than zero?
dx = 0.05;

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
    B{iD} = Data{iD} - expectedA(Data{iD});  % also catch expectedA errors here
    
    % get eigenvalues & eigenvectors 
    [V,egs] = eig(B{iD},'vector');
    [egs,ix] = sort(egs,'descend');    % sort into descending order
    V = V(:,ix);                       % ditto the eigenvectors 
 
    % get embedding dimensions
    ixRetain = find(egs > egmin); 
    M = min(4,numel(ixRetain)+1);  % how many groups? Gives huge numbers for noise, so force to tiny set of dimensions
    ixRetain = 1:M-1;
    EmbedD = V(:,ixRetain);
    
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
                    if c ~= (Mlist(iM)-1)*Treps(iT) || r ~= n
                        error('K-means sweep output is the wrong size')
                    end
                end
                if Mlist(iM) == M && Treps(iT) == max(Treps) && any(strfind(dims{iOpt},'all'))
                    Clustering{iD} = allgrps;  % use this for all next tests
                end
            end
        end
    end
    
    %% Q computation
    m = sum(sum(Data{iD}))/2;  % sum of unique weights
    try
        Q = computeQ(Clustering{iD}(:,1),B{iD},m);
    catch ME
        msg = sprintf(['\n Data' num2str(iD) ', Q compute Error: \n' ME.message ]); 
        disp(msg)
    end
    
    %% Make consensus matrix
    try
        Ccons{iD} = makeConsensusMatrix(Clustering{iD});
    catch ME
        msg = sprintf(['\n Data' num2str(iD) ', make consensus matrix Error: \n' ME.message ]); 
        disp(msg)
    end
    
    %% check convergence
    try
        [blnTrans,grpscon{iD}] = CheckConvergenceConsensus(Ccons{iD});
    catch ME
        msg = sprintf(['\n Data' num2str(iD) ', check convergence Error: \n' ME.message ]); 
        disp(msg)
    end
    
    %% embed consensus
    try
        [newD{iD},newB{iD},newM{iD},egsNorm{iD}] = EmbedConsensus(Ccons{iD});
    catch ME
        msg = sprintf(['\n Data' num2str(iD) ', embed consensus Error: \n' ME.message ]); 
        disp(msg)
    end
    x = min(egsNorm{iD})-dx:dx:max(egsNorm{iD})+dx;
    figure
    hist(egsNorm{iD},x);
    title(['Data' num2str(iD) ' normalised eigenvalues']);
   
end
