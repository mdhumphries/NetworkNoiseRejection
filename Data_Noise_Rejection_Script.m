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
fontsize = 6;

if blnViz
    % Traud Mucha Porter visualisation tools
    addpath('Traud_Mucha_Porter_CommunityVisualisation/');

    % needs MATLAB BGL Toolbox on your path - change to your local path
    % here:
    if ismac
        bglpath = genpath('/Users/mqbssmhg/Dropbox/My Toolboxes/Graph_theory/matlab_bglOSX64/');  % generate path to local BGL and all its subdirectories
    else
        bglpath = genpath('C:\Users\mqbssmhg.DS\Dropbox\My Toolboxes\Graph_theory\matlab_bgl\');
    end
    % add to current MATLAB path
    addpath(bglpath); 
end


% analysis parameters
N = 100;        % repeats of permutation
alpha = 0;  % confidence interval on estimate of maxiumum eigenvalue for null model

% WCM model options
WCMOptions.Expected = 1;
WCMOptions.NoLoops = 1;

% Poisson model options
PoissOptions.Expected = 1;
PoissOptions.NoLoops = 1;

% NodeRejection options
options.Weight = 'linear'; % 'linear' is default
options.Norm = 'L2'; % L2 is default

%% load data-file
fname = 'Lesmis.mat'; 
% % load Newman network data
load(['Networks/' fname]); 
C = 1;  % conversion factor of 1
A = full(Problem.A);

% Generate node labels for later visualisation to work
nodelabels = Problem.aux.nodename;


% % SBM generation
% A_SBM = test_noise_rejection_planted_noise(50,2,'low',0.2);
% A = A_SBM.adjacency;
% nodelabels = num2str(A_SBM.membership);
% A = round(A);

% % Star Wars: special case for pre-processing
% fname = StarWarsNetworkEp5.mat;
% load(['Networks/' fname]); 
% A = StarWars.A;
% nodelabels = StarWars.Nodes;
% nodelabels = nodelabels';

% COSYNE: special case for pre-processing
% fname = 'cosyneFinalData.mat';
% load(['Networks/' fname]); 
% A = adjMatrix;
% m = cellfun('length',cosyneData.authorHash);
% nodelabels = [];
% for i = 1:numel(cosyneData.authorHash)
%     l = numel(cosyneData.authorHash{i});
%     nodelabels = [nodelabels; cosyneData.authorHash{i} blanks(max(m) - l)];
% end

%% process adjacency matrix
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
set(gca,'Ytick',1:length(A));
set(gca,'Yticklabel',nodelabels,'Fontsize',fontsize);
% set(gca,'XTickLabelRotation',90);

%% get expected distribution of eigenvalues under null model (here, WCM)

% [Emodel,diagnostics,Vmodel] = WeightedConfigModel(A,N);

% [Emodel,diagnostics,Vmodel,ExpWCM] = WeightedConfigModel(A,N,C,WCMOptions);

[Emodel,diagnostics,Vmodel,ExpWCM] = RndPoissonConfigModel(A,N,C,PoissOptions);


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
    txt(i) = text(xynew(i,1),xynew(i,2),nodelabels(i,:),'Fontsize',fontsize);%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
end

%% Remove node labels
% for i = 1:length(A); txt(i).Color = 'none';end

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
    text(i,sorted_norms(i),nodelabels(SNIdx(i),:),'Fontsize',fontsize);%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
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
% degree_A = sum(A);

% % Z-score
% degree_A = zscore(degree_A);
% 
% % Sort and plot
% [sort_degree_A,degreeIdx_A] = sort(degree_A);
% figure
% subplot(1,2,1);
% plot(sort_degree_A)
% title('A degree distribution')
% for i = 1:length(A); 
%     text(i,sort_degree_A(i),nodelabels(degreeIdx_A(i),:),'Fontsize',fontsize);%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
% end

% Degree of nodes in modularity matrix
% degree_B = sum(B);
% 
% % Z-score
% degree_B = zscore(degree_B);
% 
% % Sort and plot
% [sort_degree_B,degreeIdx_B] = sort(degree_B);
% subplot(1,2,2);
% plot(sort_degree_B);
% title('B degree distribution')
% for i = 1:length(A); 
%     text(i,sort_degree_B(i),nodelabels(degreeIdx_B(i),:),'Fontsize',fontsize);%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
% end

%% Directly compare degree with lowD projection norm
% figure; hold all
% for i = 1:length(A);
%     plot([1,2,3],[degree_A(i),R.Difference.Norm(i),degree_B(i)],'color',[0.5,0.5,0.5])
%     text(3,degree_B(i),nodelabels(i,:),'Fontsize',fontsize);
% end
% plot(ones(length(A),1),degree_A,'.')
% plot(2*ones(length(A),1),R.Difference.Norm,'.')
% plot(3*ones(length(A),1),degree_B,'.')
% set(gca,'XTick',[1,2,3])
% set(gca,'XTickLabels',{'A degree','Eig projection','B degree'})
% xlim([0,4])

%% Scatter plot of degree vs lowD projection norm
% figure;hold all
% colormap lines
% cmap = colormap;
% 
% plot(R.Difference.Norm,degree_A,'.','color',cmap(1,:))
% plot(R.Difference.Norm,degree_B,'.','color',cmap(2,:))

% %% Scatter plot of A degree vs B degree
% clf;
% plot(degree_A,degree_B,'.','color',cmap(3,:))


%% analyse new signal matrix

% % first: remove any nodes without connections
kAsignal = sum(Asignal>0);
Connected.ixConnectedSignal = R.ixSignal(kAsignal > 1);  % more than 1 link
Connected.A = A(Connected.ixConnectedSignal,Connected.ixConnectedSignal);  % subset of original matrix

if blnViz
    % visualise Aconnected
    colors = repmat([0 0 0],n,1)+0.6;
    colors(Connected.ixConnectedSignal,:) = 0;  % black for signal, connected 
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
[Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,Ngrps,~] = allevsplitConTransitive(Connected.A);

% Louvain algorithm
[Connected.LouvCluster,Connected.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Connected.A,5,1);  % run 5 times; return 1st level of hierarchy only

%% plot sorted into group order

[H,Ix] = plotClusterMap(Connected.A,Connected.ConsCluster,[],'S');
title('Consensus clustering')
plotorder = Connected.ixConnectedSignal(Ix);

% Add node labels
numConnected = length(Connected.ixConnectedSignal);
set(gca,'Ytick',1:numConnected);
set(gca,'Yticklabel',nodelabels(plotorder,:),'Fontsize',fontsize);
% set(gca,'XTickLabelRotation',90);

for i=1:numel(Connected.LouvCluster)
    CLou = Connected.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
    [H,Ix] = plotClusterMap(Connected.A,CLou,[],'S');
    plotorder = Connected.ixConnectedSignal(Ix);
    title(['Louvain ' num2str(i)]);
    % Add node labels
    set(gca,'Ytick',1:numConnected);
    set(gca,'Yticklabel',nodelabels(plotorder,:),'Fontsize',fontsize);
    % set(gca,'XTickLabelRotation',90);
    for j = i+1:numel(Connected.LouvCluster)
        CLou2 = Connected.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
        Connected.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log(numConnected);
    end
end

%% without noise rejection
[Full.QmaxCluster,Full.Qmax,Full.ConsCluster,Full.ConsQ,Ngrps,~] = allevsplitConTransitive(A);

[H,Ix] = plotClusterMap(A,Full.ConsCluster,[],'S');
title('Consensus clustering of all')
plotorder = Ix;

% Add node labels
[srt,I] = sort(Full.ConsCluster,'ascend');
set(gca,'Ytick',1:numel(nz_e));
set(gca,'Yticklabel',nodelabels(plotorder,:),'Fontsize',fontsize);


% Louvain algorithm
[Full.LouvCluster,Full.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(A,5,1);  % run 5 times; return 1st level of hierarchy only

for i=1:numel(Full.LouvCluster)
    CLou = Full.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
    [HL,Ix] = plotClusterMap(A,CLou,[],'S');
    title(['Full Louvain ' num2str(i)]);
    plotorder = Ix;
    % Add node labels
    % [srt,I] = sort(CLou,'ascend');
    set(gca,'Ytick',1:numel(nz_e));
    set(gca,'Yticklabel',nodelabels(plotorder,:),'Fontsize',fontsize);
    for j = i+1:numel(Full.LouvCluster)
        CLou2 = Full.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
        Full.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log(numel(nz_e));
    end
end

%% save
save(['Rejected_' fname],'R','Full','Connected','A','nodelabels')
