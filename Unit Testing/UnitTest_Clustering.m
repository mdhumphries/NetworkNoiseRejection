% unit-test consensus clustering of network functions
clear all; close all
addpath ../

egmin = 1e-2; % counts as more than zero?
dx = 0.05;

%% parameter options
dims = {'all','scale'};
Treps = [0,1,10];

load UnitTestClusteringData.mat

%% test TestData 
for iD = 1:numel(TestData)
    n = size(TestData(iD).W,1); 
    % get modularity: do this with configuration model to get appropriate TestData
    P{iD} = expectedA(TestData(iD).W);
    B{iD} = TestData(iD).W - P{iD};  % also catch expectedA errors here
    
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
    Ulist = [0,1,U(iD),U(iD)+5];  % test passing of no groups, 1 group, matching groups, and more than dimensions
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
    m = sum(sum(TestData(iD).W))/2;  % sum of unique weights
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
    
    %% consensus null model
    K = 1+U(iD)-L;
    pars = {1,ones(K,1)*max(Treps),n};  % must throw error
    [blnExpected,~,~] = doUnitTest('nullmodelConsensusSweep',pars);
    pars = {L:U(iD),ones(K,1)*max(Treps),n}; 
    [blnExpected,PconS{iD},VconS{iD}] = doUnitTest('nullmodelConsensusSweep',pars);
    
    pars = {Ccons{iD}};
    [blnExpected,PconE{iD}] = doUnitTest('nullmodelConsensusExpectation',pars);

    %% embed consensus: Wishart
    pars = {Ccons{iD}};
    [blnExpected,newD{iD},newB{iD},newM{iD},egsNorm{iD}] = doUnitTest('EmbedConsensusWishart',pars);
    
    %% embed consensus: Laplacian
    pars = {Ccons{iD},U(iD)};
    [blnExpected,newDlw{iD},egslw{iD}] = doUnitTest('ProjectLaplacian',pars);
    
    %% embed consensus: using null models
    pars = {Ccons{iD},'expect'};
    [blnExpected,Dexpect{iD},Bexpect{iD},Mexpect(iD)] = doUnitTest('EmbedConsensusNull',pars);
    
    pars = {Ccons{iD},'sweep',L:U(iD),ones(K,1)*max(Treps)};
    [blnExpected,Dsweep{iD},Bsweep{iD},Msweep(iD)] = doUnitTest('EmbedConsensusNull',pars);
    
%     x = min(egsNorm{iD})-dx:dx:max(egsNorm{iD})+dx;
%     figure
%     hist(egsNorm{iD},x);
%     title(['TestData' num2str(iD) ' normalised eigenvalues']);
   
    %% entire function
    [grps{iD},Qmax(iD),grpscon{iD},Qcon(iD),ctr(iD),CLU{iD}] = ConsensusCommunityDetect(TestData(iD).W,P{iD},L,U(iD),Treps(end),dims{1});
    % same lower and upper bounds
    [grpsLU{iD},QmaxLU(iD),grpsconLU{iD},QconLU(iD),ctrLU(iD),CLU_LU{iD}] = ConsensusCommunityDetect(TestData(iD).W,P{iD},U(iD),U(iD),Treps(end),dims{1});

end

%% inspect Q, grps etc for consistency with perfect answers...

