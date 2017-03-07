%% make Network Data for Unit-testing

clear all; close all;

%% make data
% uniformly random weighted network
n = 100;
W = rand(100); % uniformly random weight matrix
W = (W+W')/2;    % undirected 
W(eye(n)==1) = 0; % remove self-connections
Data{1} = W;

% block diagoonal
Wblock = ones(n/2); 
W = zeros(n);
W(1:n/2,1:n/2) = Wblock;
W(n/2+1:end,n/2+1:end) = Wblock;
W(eye(n)==1) = 0; % remove self-connections
Data{2} = W;

% multiple groups
Wblock = ones(n/4); 
W = zeros(n);
for iW = 1:4
    ix = (1:n/4)+n/4*(iW-1);
    W(ix,ix) = Wblock;
end
W(eye(n)==1) = 0; % remove self-connections
Data{3} = W;

save UnitTestClusteringData Data