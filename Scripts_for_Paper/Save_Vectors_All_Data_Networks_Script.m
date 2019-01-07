%% script to save all eigenvectors from the null models of each data network
% 
% Mark Humphries 31/07/2018

clearvars;

storepath = 'C:\Users\lpzmdh\Dropbox\Analyses\Networks\DataNets_Null_EigVectors\';


networks = {'StarWarsNetworkEp1','StarWarsNetworkEp2','StarWarsNetworkEp3','StarWarsNetworkEp4','StarWarsNetworkEp5','StarWarsNetworkEp6',...
            'Allen_Gene_Leaf','Lesmis','adjnoun','polblogs','dolphins','cosyneFinalData','celegansneural','power','cElegAdjMatAllSynapUndirected'};

nNet = numel(networks);
% sanity check: run through list, and check all can be loaded without error
for iN = 1:nNet
    load(['../Networks/' networks{iN}]);
end

addpath('../Network_Spectra_Functions/')
addpath('../Network_Analysis_Functions/')

%% get null model eigenvectors for each network, and store
for iN = 1:nNet

    fname = networks{iN}
    % load rejected network to get processed original weight matrix in
    % Data.A
    % and to get all analysis parameters...
    
    load(['../Results/Rejected_' fname])
    
    % run rejection
    switch pars.Model
    case 'Poiss'
        [Emodel,diagnostics,Vmodel] = poissonSparseWCM(Data.A,pars.N,pars.C,optionsModel);
    case 'Link'
        [Emodel,diagnostics,Vmodel] = linkSparseWCM(Data.A,pars.N,pars.C,optionsModel);
    otherwise
        error('Unrecognised null model specified')
    end
    
    % save eigenvectors of null model (may as well save eigenvalues too...)
    save([storepath '/NullModel_Eigenspectrum_' fname],'Emodel','Vmodel','pars','-v7.3')
end