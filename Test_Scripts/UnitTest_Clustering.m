% unit-test consensus clustering of network functions
clear all; close all
addpath ../

fcns = {'kmeansSweep','computeQ','makeConsensusMatrix','CheckConvergenceConsensus','EmbedConsensus'};

egmin = 1e-2; % counts as more than zero?
dx = 0.05;

%% parameter options
dims = {'all','scale'};
Treps = [0,1,10];

load UnitTestClusteringData

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
    L = 2;
    U = min(4,numel(ixRetain)+1);  % how many groups? Gives huge numbers for noise, so force to tiny set of dimensions
    ixRetain = 1:U-1;
    EmbedD = V(:,ixRetain);
    
    %% test k-means sweep
    Ulist = [0,1,U,U+5];  % test passing of no groups, 1 group, matching groups, and more than dimensions
    ErrCtr = 1;
    for iM = 1:numel(Ulist)
        for iT = 1:numel(Treps)
            for iOpt = 1:numel(dims)
                pars = {EmbedD,L,Ulist(iM),Treps(iT),dims{iOpt}};
                [blnExpected,allgrps] = doUnitTest(fcns{1},pars);
                if blnExpected
                    allgrps = [];
                end
%                 try
%                     allgrps = kmeansSweep(EmbedD,L,Ulist(iM),Treps(iT),dims{iOpt});
%                 catch ME
%                     if any(strfind(ME.identifier,'kmeansSweep'))
%                         msg = sprintf(['\n Data' num2str(iD) ', K-means Error: \n' ME.message ' \n Thrown by: M=' num2str(Ulist(iM)) ', Treps = ' num2str(Treps(iT)) ', dims = ' dims{iOpt}]); 
%                         disp(msg)
%                         allgrps = [];
%                         ErrCtr = ErrCtr+1;
%                     else
%                        rethrow(ME); 
%                     end
%                 end
                
                
                % sanity checks
                if ~isempty(allgrps)
                    [r,c] = size(allgrps);
                    if c ~= (Ulist(iM)-1)*Treps(iT) || r ~= n
                        error('K-means sweep output is the wrong size')
                    end
                end
                if Ulist(iM) == U && Treps(iT) == max(Treps) && any(strfind(dims{iOpt},'all'))
                    Clustering{iD} = allgrps;  % use this for all next tests
                end
            end
        end
    end
    
    %% Q computation
    m = sum(sum(Data{iD}))/2;  % sum of unique weights
    pars = {Clustering{iD}(:,1),B{iD},m};
    [blnExpected,Q] = doUnitTest(fcns{1},pars);
   
    
    %% Make consensus matrix
    try
        Ccons{iD} = makeConsensusMatrix(Clustering{iD});
    catch ME
        msg = sprintf(['\n Data' num2str(iD) ', make consensus matrix Error: \n' ME.message ]); 
        disp(msg)
    end
    
    %% check convergence
    pars = Ccons(iD);
    [blnExpected,blnTrans,grpscon{iD}] = doUnitTest(fcns{2},pars);
    
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
   
    %% entire function
    [grps{iD},Qmax(iD),grpscon{iD},Qcon(iD),ctr(iD),CLU{iD}] = ConsensusCommunityDetect(Data{iD},B{iD},U,Treps(end),dims{1});
end
