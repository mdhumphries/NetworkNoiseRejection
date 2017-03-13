% make original trilogy network
clear all; close all;

Eps = [4,5,6];

% load each episode
for iE = 1:numel(Eps)
    load(['StarWarsNetworkEp' num2str(Eps(iE))]);
    Net(iE).A = StarWars.A;
    Net(iE).nodelabels = StarWars.Nodes;
end

%% merge: master node-labels list
nodelabels = Net(1).nodelabels;
for iE = 2:numel(Eps)
    for iN = 1:numel(Net(iE).nodelabels)
        L = ismember(nodelabels,Net(iE).nodelabels{iN});  % is this node already in the list?
        if ~any(L)
            nodelabels(numel(nodelabels)+1) = Net(iE).nodelabels(iN);
        end
    end
end

%% merge into one network
StarWars = [];
StarWars.A = zeros(numel(nodelabels));
StarWars.Nodes = nodelabels;

for iN = 1:numel(nodelabels)
    for iM = iN+1:numel(nodelabels)
        for iE = 1:numel(Eps)
            L1 = ismember(Net(iE).nodelabels,nodelabels(iN)); % find first node label in current network
            L2 = ismember(Net(iE).nodelabels,nodelabels(iM)); % find second node label in current network
            if any(L1) && any(L2)
                StarWars.A(iN,iM) = StarWars.A(iN,iM) + Net(iE).A(L1,L2);  % if both there, add weight
            end
        end
    end
end

StarWars.A = StarWars.A + StarWars.A'; 

save StarWarsOriginalTrilogy StarWars