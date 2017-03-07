% unit-test consensus clustering of network functions
clear all; close all
addpath ../

egmin = 1e-2; % counts as more than zero?
dx = 0.05;

%% parameter options
dims = {'all','scale'};
Treps = [0,1,10];

load UnitTestClusteringData

%% test data 
for iD = 1:numel(Data)
    n = size(Data{iD},1); 
    % get modularity: do this with configuration model to get appropriate data
    P{iD} = expectedA(Data{iD});
    B{iD} = Data{iD} - P{iD};  % also catch expectedA errors here
    
    % get eigenvalues & eigenvectors 
    [V,egs] = eig(B{iD},'vector');
    [egs,ix] = sort(egs,'descend');    % sort into descending order
    V = V(:,ix);                       % ditto the eigenvectors 
    
    % get embedding dimensions
    ixRetain = find(egs > egmin); 
    L = 2;
    U(iD) = min(4,numel(ixRetain)+1);  % how many groups? Gives huge numbers for noise, so force to tiny set of dimensions
    ixRetain = 1:U(iD)-1;
    EmbedD = V(:,ixRetain);
    
    %% test k-means sweep
    Ulist = [0,1,U(iD),U+5];  % test passing of no groups, 1 group, matching groups, and more than dimensions
    ErrCtr = 1;
    for iM = 1:numel(Ulist)
        for iT = 1:numel(Treps)
            for iOpt = 1:numel(dims)
                pars = {EmbedD,L,Ulist(iM),Treps(iT),dims{iOpt}};
                [blnExpected,allgrps] = doUnitTest('kmeansSweep',pars);
             
                % sanity checks
                if ~isempty(allgrps)
                    [r,c] = size(allgrps);
                    if c ~= (Ulist(iM)-1)*Treps(iT) || r ~= n
                        error('K-means sweep output is the wrong size')
                    end
                end
                if Ulist(iM) == U(iD) && Treps(iT) == max(Treps) && any(strfind(dims{iOpt},'all'))
                    Clustering{iD} = allgrps;  % use this for all next tests
                end
            end
        end
    end

    %% Q computation
    m = sum(sum(Data{iD}))/2;  % sum of unique weights
    pars = {Clustering{iD}(:,1),B{iD},m};
    [blnExpected,Q] = doUnitTest('computeQ',pars);
   
    pars = {Clustering{iD}(:,1),B{iD}};
    [blnExpected,Q] = doUnitTest('computeQ',pars);
   
    %% Make consensus matrix
    pars = {Clustering{iD}};
    [blnExpected,Ccons{iD}] = doUnitTest('makeConsensusMatrix',pars);
    
    %% check convergence
    pars = {Ccons{iD}};
    [blnExpected,blnTrans,grpscon{iD}] = doUnitTest('CheckConvergenceConsensus',pars);

    
    %% embed consensus
    pars = {Ccons{iD}};
    [blnExpected,newD{iD},newB{iD},newM{iD},egsNorm{iD}] = doUnitTest('EmbedConsensus',pars);

    x = min(egsNorm{iD})-dx:dx:max(egsNorm{iD})+dx;
    figure
    hist(egsNorm{iD},x);
    title(['Data' num2str(iD) ' normalised eigenvalues']);
   
    %% entire function
    [grps{iD},Qmax(iD),grpscon{iD},Qcon(iD),ctr(iD),CLU{iD}] = ConsensusCommunityDetect(Data{iD},P{iD},U(iD),Treps(end),dims{1});

end

%% inspect Q, grps etc for consistency with perfect answers...

