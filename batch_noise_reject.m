% Script to run through all the networks in the example folder, reject the
% noise nodes and save an array of signal node labels for later analysis.

network_list = {'adjnoun';'celegansneural';'dolphins';'karate';...
    'lesmis';'polblogs';'power';'cosyneFinalData';'cond-mat';'cond-mat-2005';};

for i = 7:length(network_list)
    display(['Processing ',network_list{i},'.....'])
    
    load(['Networks/',network_list{i},'.mat']);
    
    A = full(Problem.A);
    
    reject_the_noise(A,network_list{i});
end

%% Process Star Wars network separately due to different file format
load('Networks/StarWarsNetworkAll.mat')

A = StarWars.A;

A = full(A);

reject_the_noise(A,'StarWars');