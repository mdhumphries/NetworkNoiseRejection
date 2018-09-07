%% example weight distributions

clearvars; close all

addpath ../Helper_Functions/
addpath ../Network_Spectra_Functions/
addpath ../Network_Analysis_Functions/
addpath ../ZhangNewman2015/

load '../Networks/lesmis.mat';
A = full(Problem.A);

% number of trials
noTrials = 30;

% max number of groups to test
nGroups = 20;

for count = 1:noTrials
    [bestPartition,maxQPartition,Results.Q(count,:),Results.ixBest(count),Results.ixQ(count),Results.nGroups(count,:)] =...
            multiwaySpectCommDet(A,nGroups);

     Results.bestNoOfGroups(count) = max(bestPartition);
     Results.QmaxNoOfGroups(count) = max(maxQPartition);
end

figure
subplot(211),
hist(Results.bestNoOfGroups, 1:max(Results.bestNoOfGroups))
ylabel('frequency')
xlabel('No. of groups detected: knee')
subplot(212),
hist(Results.QmaxNoOfGroups, 1:max(Results.QmaxNoOfGroups))
ylabel('frequency')
xlabel('No. of groups detected: Qmax')

%% plateaus 
figure
%for iN = 1:noTrials
yyaxis left
plot(1:nGroups,Results.nGroups(:,1:end),'.-','Linewidth',0.5)
%end
xlabel('Number of groups tested')
ylabel('Number of groups found')

yyaxis right
plot(1:nGroups,Results.Q(:,1:end),'.-','Linewidth',0.5)
%end
ylabel('Q')

save ../Results/MultiwayExample Results noTrials nGroups
