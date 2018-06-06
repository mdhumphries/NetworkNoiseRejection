%% script to visualise synthetic model results
% Mark Humphries 8/8/2017

clear all; close all;

addpath('../Helper_Functions/')

fname = 'P_rejection_SyntheticEqual_NoNoise_20180606T123256';
blnCluster = 0;  % if done clustering, set this to 1

% load results and set maps
load(['../Results/' fname]);  

nP = numel(Model.P_of_within);    % within communities P
nBatch = size(Results.Time,2);   % number of model repeats

% plotting parameters
I = 0.99; % 99% CI

% colormaps
perf_cmap = brewermap(10,'Blues');
VI_cmap = brewermap(10,'Reds');
colors.line = [0.8 0.7 0.7];
colors.error = [0.8 0.7 0.7];
colors.WCM = [0.8 0.5 0.5];
colors.Config= [0.5 0.5 0.8];
colors.WCMpos = [0.9 0.7 0.7];
colors.Configpos= [0.7 0.7 0.9];
colors.truth = [0.7 0.7 0.7];

M = 5;

Pdiff = Model.P_of_within - Model.P.between;

%% processing time
meanTime = mean(Results.Time,2);
CITime = CIfromSEM(std(Results.Time,[],2),ones(nP,1).*nBatch,I);

figure
plot(Pdiff,meanTime,'o-','Color',colors.line,'MarkerSize',M,'MarkerFaceColor',colors.line)
line([Pdiff; Pdiff],[(meanTime-CITime)';(meanTime+CITime)'],'Color',colors.error)
set(gca,'XTick',Pdiff,'XTickLabel',Pdiff)
ylabel('Time (s)')
xlabel('P(within) - P(between)')
box off
%% recovered proportions of networks with groups

figure
plot(Pdiff,Results.ProportionModular.SpectraWCM ./nBatch,'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); hold on
plot(Pdiff,Results.ProportionModular.SpectraConfig ./nBatch,'o-','MarkerSize',M,'Color',colors.Config,'MarkerFaceColor',colors.Config)
plot(Pdiff,Results.ProportionModular.PosEigWCM ./nBatch,'o-','MarkerSize',M,'Color',colors.WCMpos,'MarkerFaceColor',colors.WCMpos);
set(gca,'XTick',Pdiff,'XTickLabel',Pdiff)
xlabel('P(within) - P(between)')
ylabel('Networks with modules (%)')
box off

%% number of groups
mWCM = mean(Results.SpectraWCMGroups,2);
mPosWCM = mean(Results.PosEigWCMGroups,2);
[bndsWCM.L,bndsWCM.U] = bounds(Results.SpectraWCMGroups,2);
[bndsPosWCM.L,bndsPosWCM.U] = bounds(Results.PosEigWCMGroups,2);

figure
line([0 max(Pdiff)],[numel(Model.N), numel(Model.N)],'Color',colors.truth); hold on
plot(Pdiff,mWCM,'o-','MarkerSize',M,'Color',colors.WCM,'MarkerFaceColor',colors.WCM); hold on
line([Pdiff; Pdiff],[bndsWCM.L';bndsWCM.U'],'Color',colors.WCM)

plot(Pdiff,mPosWCM,'o-','MarkerSize',M,'Color',colors.WCMpos,'MarkerFaceColor',colors.WCMpos);
line([Pdiff; Pdiff],[bndsPosWCM.L';bndsPosWCM.U'],'Color',colors.WCMpos)

set(gca,'XTick',Pdiff,'XTickLabel',Pdiff)
xlabel('P(within) - P(between)')
ylabel('Number of modules')
box off

%% VI of clustering
if blnCluster
    % load clustering
    load(['../Results/Clustering' fname]);  % Main table is 'core_perf'

end