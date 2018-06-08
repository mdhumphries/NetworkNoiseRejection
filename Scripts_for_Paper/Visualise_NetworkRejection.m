%% script to visualise Full and Signal networks
% Visualisations of network layouts need:
% (1) Traud-Mucha-Porter toolbox (included in this GitHub repo)
% (2) MATLAB BGL Toolbox:
%           Win32, Win64, Mac32, Linux: https://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/ 
%           Mac64: http://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/old/matlab_bgl_4.0_osx64.zip
%
% Mark Humphries, Mat Evans 28/2/2017

clear all; close all
addpath('../Helper_Functions/')

fname = 'StarWarsNetworkEp6'; 
% fname = 'LesMis'; 
% fname = 'cosyneFinalData';
blnVizNet = 0;  % network visualisation - if MATLAB BGL installed, appropriate for platform:
blnExport = 1;

blnLabel = 1;  % add name labels to matrix plots
fontsize = 6;
Nodes = 5; % marker size for nodes
Links = 1; % linewidth for links

%% load data
load(['../Results/Rejected_' fname],'Data','Rejection')

%% ready code for visualising networks
if blnVizNet
    % Traud Mucha Porter visualisation tools
    addpath('../Traud_Mucha_Porter_CommunityVisualisation/');

    % needs MATLAB BGL Toolbox on your path - *change to your local path here*:
    if ismac
        bglpath = genpath('/Users/mqbssmhg/Dropbox/My Toolboxes/Graph_theory/matlab_bglOSX64/');  % generate path to local BGL and all its subdirectories
    else
        bglpath = genpath('C:\Users\lpzmdh\Dropbox\My Toolboxes\Graph_theory\matlab_bgl\');
    end
    % add to current MATLAB path
    addpath(bglpath); 
end


%% Image plot of A, with labels
k = sum(Data.A);
[~,iK] = sort(k,'descend');   % plot sorted into degree order
cmap = brewermap(10,'Greys');
figure; 
imagesc(Data.A(iK,iK)); colormap(cmap);
if blnLabel 
    set(gca,'Ytick',1:length(Data.ixRetain)); 
    set(gca,'Yticklabel',Data.nodelabels(iK,:),'Fontsize',fontsize); 
end
if blnExport exportPPTfig(gcf,[fname 'W'],[10 15 25 20]); end

%% compare eigenvalues of data and model
B = Data.A - Data.ExpA;  % modularity matrix using chosen null model
Edata = eig(B);
[fdata,xdata] = ecdf(Edata);
[fmodel,xmodel] = ecdf(Data.Emodel(:));
figure
stairs(xdata,fdata,'r'); hold on
stairs(xmodel,fmodel,'k')
ylabel('P')
xlabel('Eigenvalues')

%% compare data eigenvalues to distribution of maximum eigenvalues
if isfield(Data,'NEigEst')
maxModelE = max(Data.Emodel);   % all predicted maximum values
minModelE = min(Data.Emodel);   % all predicted minimum values

[F_maxmodelE,X_maxmodelE] = ksdensity(maxModelE);
%CImax = CIfromSEM(std(maxModelE),numel(maxModelE),[0.95,0.99]);

[F_minmodelE,X_minmodelE] = ksdensity(minModelE);
% ExpMin = mean(minModelE);
%CImin = CIfromSEM(std(minModelE),numel(minModelE),[0.95,0.99]);


figure
% density histograms of min and max
plot(X_maxmodelE,F_maxmodelE,'color',[0.7 0.7 0.7]); hold on
plot(X_minmodelE,F_minmodelE,'color',[0.7 0.7 0.7]); hold on
ylim = get(gca,'YLim');

% compute PIs
line([Data.EigEst(1), Data.EigEst(1)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1); % lower CI (or only CI, if just computing mean)
line([Data.EigEst(2), Data.EigEst(2)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1); % upper CI (or 0, if just computing mean)

line([Data.NEigEst(1), Data.NEigEst(1)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1); % lower CI
line([Data.NEigEst(2), Data.NEigEst(2)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1); % upper CI


% plot CIs and means
% line([Data.EigEst(1), Data.EigEst(1)],ylim,'Color',[0.8 0.6 0.6],'Linewidth',2); % mean of max Eigs
% % 95% CIs
% line([Data.EigEst(1)-CImax(1), Data.EigEst(1)-CImax(1)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1); % lower CI
% line([Data.EigEst(1)+CImax(1), Data.EigEst(1)+CImax(1)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1); % upper CI
% % 99% CIs
% line([Data.EigEst(1)-CImax(2), Data.EigEst(1)-CImax(2)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1,'Linestyle',':'); % lower CI
% line([Data.EigEst(1)+CImax(2), Data.EigEst(1)+CImax(2)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1,'Linestyle',':'); % upper CI

% and the same for the minimum
% line([ExpMin, ExpMin],ylim,'Color',[0.8 0.6 0.6],'Linewidth',2); % mean of min Eigs
% % 95% CIs
% line([ExpMin - CImin(1), ExpMin - CImin(1)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1); % lower CI
% line([ExpMin + CImin(1), ExpMin + CImin(1)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1); % upper CI
% % 99% CIs
% line([ExpMin - CImin(2), ExpMin - CImin(2)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1,'Linestyle',':'); % lower CI
% line([ExpMin + CImin(2), ExpMin + CImin(2)],ylim,'Color',[0.6 0.6 0.8],'Linewidth',1,'Linestyle',':'); % upper CI



% data
plot(Edata,zeros(numel(Edata),1)+0.01,'x')

title('Max/min model eigenvalues vs data')
ylabel('P(max model eigenvalue)')
xlabel('Eigenvalue')

end

%% just plot eigenvalues and limits

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2]); 
line([Data.EigEst(1), Data.EigEst(1)],[0 0.02],'Color',[0.6 0.6 0.8],'Linewidth',1); hold on % lower CI
line([0,0],[0 0.02],'Color',[0.8 0.8 0],'Linewidth',1); % modularity matrix
plot(Edata,zeros(numel(Edata),1)+0.01,'kx')
lim = max(max(abs(Edata)),Data.EigEst(1)) + 2;
set(gca,'XLim',[-lim lim])
% ylabel('P(max model eigenvalue)')
xlabel('Eigenvalues')
if blnExport exportPPTfig(gcf,[fname 'Eigs_Limits'],[10 15 4 2]); end

%% visualise signal and noise parts
if blnVizNet
    n = size(Data.A,1);
    xynew = fruchterman_reingold_force_directed_layout(sparse(Data.A));
    syms = repmat('o',n,1);
    colors = repmat([0 0 0],n,1); 
    figure
    graphplot2D(xynew,Data.A,1,colors,syms,Nodes);
    axis off
    allh = get(gca,'Children');
    set(allh,'MarkerEdgeColor',[0 0 0])
    set(allh,'LineWidth',Links)
    %title('Full')
    if blnExport exportPPTfig(gcf,[fname 'FullNetwork'],[10 15 25 25]); end
    
    % now shade out rejection
    colors(Rejection.ixNoise,:) = colors(Rejection.ixNoise,:) + 0.6;  % gray for noise 
    figure
    graphplot2D(xynew,Data.A,1,colors,syms,Nodes);
    axis off
    allh = get(gca,'Children');
    set(allh,'MarkerEdgeColor',[0 0 0])
    set(allh,'LineWidth',Links)
    title('Signal/Noise split')
    % Add node labels to graph
    if blnLabel
        for i = 1:length(Data.ixRetain); 
            txt(i) = text(xynew(i,1),xynew(i,2),Data.nodelabels(i,:),'Fontsize',fontsize);%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
        end
    end
end



%% Plot lowD projection of each node, with labels, sorted by magnitude
if Data.Dn > 0
    
    figure
    [sorted_norms,SNIdx] = sort(Rejection.Difference.Norm); % SNIdx = Sorted Norm Index

    stem(1:numel(Rejection.ixNoise),sorted_norms(1:numel(Rejection.ixNoise)),'Color',[0.7 0.7 0.7]);
    hold all
    stem(numel(Rejection.ixNoise)+1:numel(Rejection.Difference.Norm),sorted_norms(numel(Rejection.ixNoise)+1:end),'Color',[0.8 0.5 0.55])
    xlabel('Nodes')
    ylabel('Projection (normalised to null model)')

    % EDIT HERE: text labels 90 rotated: up for y>0; down for y<0
    for i = 1:length(Data.A); 
        text(i,sorted_norms(i),Data.nodelabels(SNIdx(i),:),'Fontsize',fontsize);%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
    end
    
    
    figure
    [sorted_norms,SNIdx] = sort(Rejection.NegDiff.Norm); % SNIdx = Sorted Norm Index

    stem(1:numel(Rejection.ixNegative),sorted_norms(1:numel(Rejection.ixNegative)));
    hold all
    stem(numel(Rejection.ixNegative)+1:numel(Rejection.NegDiff.Norm),sorted_norms(numel(Rejection.ixNegative)+1:end),'Color',[0.7 0.7 0.7]);
    xlabel('Nodes')
    ylabel('Projection (normalised to null model)')

    % EDIT HERE: text labels 90 rotated: up for y>0; down for y<0
    for i = 1:length(Data.A); 
        text(i,sorted_norms(i),Data.nodelabels(SNIdx(i),:),'Fontsize',fontsize);%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
    end

else
    warning('No retained dimensions, so no node rejection either')
end

%% visualise connected-signal network
if blnVizNet
    % visualise Aconnected
    colors = repmat([0 0 0],n,1)+0.6;
    colors(Data.ixSignal_Final,:) = 0;  % black for signal, connected 
    figure
    graphplot2D(xynew,Data.A,1,colors,syms,Nodes);
    axis off
    allh = get(gca,'Children');
    set(allh,'MarkerEdgeColor',[0 0 0])
    set(allh,'LineWidth',Links)
    % title('Connected Signal/Noise split')
    if blnExport exportPPTfig(gcf,[fname 'FinalSignalNetwork'],[10 15 8 8]); end
    
    % Add node labels to graph
    for i = Data.ixSignal_Final; 
        txt(i) = text(xynew(i,1),xynew(i,2),Data.nodelabels(i,:),'Fontsize',fontsize);%,'BackgroundColor',[0.9,0.9,0.9],'alpha',0.5);
    end
end




