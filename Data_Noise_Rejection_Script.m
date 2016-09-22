%% template script for applying complete work-flow to one data network, using one choice of null model
% data network: correlations between firing in Aplysia recording
% null model: weighted configuration model
%
% Visualisations need:
% (1) Traud-Mucha-Porter toolbox (included in GitHub)
% (2) MATLAB BGL Toolbox:
%           Win32, Win64, Mac32, Linux: https://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/ 
%           Mac64: http://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/old/matlab_bgl_4.0_osx64.zip

clear all; close all
blnViz = 1;  % if MATLAB BGL installed, appropriate for platform:

if blnViz
    % Traud Mucha Porter visualisation tools
    addpath('Traud_Mucha_Porter_CommunityVisualisation/');

    % needs MATLAB BGL Toolbox on your path - change to your local path
    % here:
    % bglpath = genpath('/Users/mqbssmhg/Dropbox/My Toolboxes/Graph_theory/matlab_bglOSX64/');  % generate path to local BGL and all its subdirectories
    bglpath = genpath('C:\Users\mqbssmhg.DS\Dropbox\My Toolboxes\Graph_theory\matlab_bgl\');
    
    % add to current MATLAB path
    addpath(bglpath); 
end


% analysis parameters
N = 100;        % repeats of permutation
alpha = 0;  % confidence interval on estimate of maxiumum eigenvalue for null model

% WCM model options
WCMOptions.Expected = 1;
WCMOptions.NoLoops = 1;

% NodeRejection options
options.Weight = 'linear'; % 'linear' is default
options.Norm = 'L2'; % L2 is default

% % load Newman network data
load('Networks/Lesmis.mat');
% % load('Networks/dolphins.mat');
% load('Networks/polblogs.mat');
A = full(Problem.A);

% Generate node labels for later visualisation to work
nodelabels = Problem.aux.nodename;


% % SBM generation
% A_SBM = test_noise_rejection_planted_noise(50,2,'low',0.2);
% A = A_SBM.adjacency;
% nodelabels = num2str(A_SBM.membership);
% A = round(A);

% % Star Wars
% load('Networks/StarWarsNetworkAll.mat')
% A = StarWars.A;
% nodelabels = StarWars.Nodes;
% nodelabels = nodelabels';

% Ensure no self loops
A(find(eye(length(A)))) = 0;

% Strip out zero elements
nz_e = find(sum(A)); % nonzero_elements
A = A(nz_e,nz_e);
nodelabels = nodelabels(nz_e,:);

% Ensure A is symetric
A = (A + A')/2;

% Image plot of A, with labels
imagesc(A);
set(gca,'Xtick',1:length(A));
set(gca,'Xticklabel',nodelabels);
set(gca,'XTickLabelRotation',90);

% get expected distribution of eigenvalues under null model (here, WCM)

% [Emodel,diagnostics,Vmodel] = WeightedConfigModel(A,N);

[Emodel,diagnostics,Vmodel,ExpWCM] = WeightedConfigModel(A,N,1,WCMOptions);


%% decompose nodes into signal and noise
% B = A - expectedA(A);  % modularity matrix using chosen null model
B = A - ExpWCM;  % modularity matrix using chosen null model

% compare data and model
Edata = eig(B);
[fdata,xdata] = ecdf(Edata);
[fmodel,xmodel] = ecdf(Emodel(:));
figure
stairs(xdata,fdata,'r'); hold on
stairs(xmodel,fmodel,'k')
xlabel('P')
ylabel('Eigenvalues')

% find low-dimensional projection
[Dspace,Ix,Dn,EigEst] = LowDSpace(B,Emodel,alpha); % to just obtain low-dimensional projection

% node rejection within low-dimensional projection
R = NodeRejection(B,Emodel,alpha,Vmodel,options); % N.B. also calls function to find projections

% new signal matrix
Asignal = A(R.ixSignal,R.ixSignal);


%% visualise signal and noise parts
if blnViz
    n = size(A,1);
    xynew = fruchterman_reingold_force_directed_layout(sparse(A));
    syms = repmat('o',n,1);

    colors = repmat([0 0 0],n,1);
    colors(R.ixNoise,:) = colors(R.ixNoise,:) + 0.6;  % gray for noise 
    figure
    graphplot2D(xynew,A,10,colors,syms,10);
    axis off
    allh = get(gca,'Children');
    set(allh,'MarkerEdgeColor',[0 0 0])
    set(allh,'LineWidth',0.2)
    title('Signal/Noise split')
end

%% Add node labels to graph
for i = 1:length(A); 
    txt(i) = text(xynew(i,1),xynew(i,2),nodelabels(i,:));%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
end

%% Remove node labels
for i = 1:length(A); txt(i).Color = 'none';end

%% Plot lowD projection of each node, with labels, sorted by magnitude
figure
[sorted_norms,SNIdx] = sort(R.Difference.Norm); % SNIdx = Sorted Norm Index

stem(1:numel(R.ixNoise),sorted_norms(1:numel(R.ixNoise)));
hold all
stem(numel(R.ixNoise)+1:numel(R.Difference.Norm),sorted_norms(numel(R.ixNoise)+1:end))
xlabel('Nodes')
ylabel('Projection (normalised to null model)')

% EDIT HERE: text labels 90 rotated: up for y>0; down for y<0
for i = 1:length(A); 
    text(i,sorted_norms(i),nodelabels(SNIdx(i),:));%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
end



% Unsorted version
% figure
% stem(R.ixSignal,R.Difference.Norm(R.ixSignal))
% hold all
% stem(R.ixNoise,R.Difference.Norm(R.ixNoise))
% 
% for i = 1:length(A); 
%     text(i,R.Difference.Norm(i),Problem.aux.nodename(i,:));%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
% end

%% Compare to node degree - no correlation
% Degree of original adjacency matrix
degree_A = sum(A);

% Z-score
degree_A = zscore(degree_A);

% Sort and plot
[sort_degree_A,degreeIdx_A] = sort(degree_A);
figure
subplot(1,2,1);
plot(sort_degree_A)
title('A degree distribution')
for i = 1:length(A); 
    text(i,sort_degree_A(i),nodelabels(degreeIdx_A(i),:));%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
end

% Degree of nodes in modularity matrix
degree_B = sum(B);

% Z-score
degree_B = zscore(degree_B);

% Sort and plot
[sort_degree_B,degreeIdx_B] = sort(degree_B);
subplot(1,2,2);
plot(sort_degree_B);
title('B degree distribution')
for i = 1:length(A); 
    text(i,sort_degree_B(i),nodelabels(degreeIdx_B(i),:));%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
end

%% Directly compare degree with lowD projection norm
figure; hold all
for i = 1:length(A);
    plot([1,2,3],[degree_A(i),R.Difference.Norm(i),degree_B(i)],'color',[0.5,0.5,0.5])
    text(3,degree_B(i),nodelabels(i,:));
end
plot(ones(length(A),1),degree_A,'.')
plot(2*ones(length(A),1),R.Difference.Norm,'.')
plot(3*ones(length(A),1),degree_B,'.')
set(gca,'XTick',[1,2,3])
set(gca,'XTickLabels',{'A degree','Eig projection','B degree'})
xlim([0,4])

%% Scatter plot of degree vs lowD projection norm
figure;hold all
colormap lines
cmap = colormap;

plot(R.Difference.Norm,degree_A,'.','color',cmap(1,:))
plot(R.Difference.Norm,degree_B,'.','color',cmap(2,:))

% %% Scatter plot of A degree vs B degree
% clf;
% plot(degree_A,degree_B,'.','color',cmap(3,:))


%% analyse new signal matrix

% % first: remove any nodes without connections
kAsignal = sum(Asignal>0);
ixConnectedSignal = R.ixSignal(kAsignal > 1);  % more than 1 link
Aconnected = A(ixConnectedSignal,ixConnectedSignal);  % subset of original matrix

if blnViz
    % visualise Aconnected
    colors = repmat([0 0 0],n,1)+0.6;
    colors(ixConnectedSignal,:) = 0;  % black for signal, connected 
    figure
    graphplot2D(xynew,A,10,colors,syms,10);
    axis off
    allh = get(gca,'Children');
    set(allh,'MarkerEdgeColor',[0 0 0])
    set(allh,'LineWidth',0.5)
    title('Connected Signal/Noise split')
end

% consensus modularity
% [C,Qmax,Ccon,Qc,Ngrps,Q] = allevsplitConTransitive(Asignal);
[C,Qmax,Ccon,Qc,Ngrps,~] = allevsplitConTransitive(Aconnected);
% [C2,Qmax2,Ccon2,Qc2,Ngrps,~] = allevsplitConTransitive(Aconnected);

% [Cfull,Qmaxfull,Cconfull,Qcfull,Nfull,~] = allevsplitConTransitive(A);

% Louvain algorithm
[allC,allQ,allCn,allIters] = LouvainCommunityUDnondeterm(Aconnected,5,1);  % run 5 times; return 1st level of hierarchy only

%% plot sorted into group order

[H,Ix] = plotClusterMap(Aconnected,Ccon,[],'S');
title('Consensus clustering')
plotorder = ixConnectedSignal(Ix);

% Add node labels
numConnected = length(Aconnected);
[srt,I] = sort(Ccon,'ascend');
set(gca,'Ytick',1:numConnected);
set(gca,'Yticklabel',nodelabels(plotorder,:));
% set(gca,'XTickLabelRotation',90);

for i=1:numel(allC)
    CLou = allC{i}{1};  % Repeat#, Level of Hierarchy
    HL = plotClusterMap(Aconnected,CLou);
    title(['Louvain ' num2str(i)]);
    % Add node labels
    [srt,I] = sort(CLou,'ascend');
    set(gca,'Ytick',1:numConnected);
    set(gca,'Yticklabel',nodelabels(plotorder,:));
    % set(gca,'XTickLabelRotation',90);
end