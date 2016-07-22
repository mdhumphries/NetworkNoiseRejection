% test deconstruction....
clear all; close all

addpath('../')

%% test on random matrices
T = 1000; % time-steps
N = 100; % number of time-series
r = rand(T,N); % random time-series
rN = (r-mean(r(:)))./std(r(:));  % z-scored
Cxy = corrcoef(r);
% Cxy = corrcoef(rN);

[Cs,Cn,ixRetain,sigR,noiseR] = deconstructCxy(Cxy,T); % ixRetain & sigR is empty, as expected

figure
subplot(3,3,1), imagesc(Cxy); title('Cxy')
subplot(3,3,2), imagesc(Cs); title('Signal'); colorbar
subplot(3,3,3), imagesc(Cn); title('Noise')

subplot(3,3,7), imagesc(Cxy(ixRetain,ixRetain)); title('Cxy retained')
subplot(3,3,8), imagesc(sigR); title('Signal, retained'); colorbar
subplot(3,3,9), imagesc(noiseR); title('Noise, retained')

% also get eigenvalues
[V,D] = eig(Cxy);
egs =diag(D);

% upper limit of random matrix distribution - Marcenko-Pastur distribution
% for IID random matrices, derived from time-series
MP_lambda_pos =(1+sqrt(N/T)).^2;
MP_lambda_neg =(1-sqrt(N/T)).^2;

subplot(3,3,[4:6])
hist(egs,30)
line([MP_lambda_neg MP_lambda_neg],[0 10],'Color',[1 0 0])  % lower M-P bound
line([MP_lambda_pos MP_lambda_pos],[0 10],'Color',[1 0 0])  % upper M-P bound
xlabel('Eigenvalues')
ylabel('Frequency')
title('Eigenvalue distribution of Cxy of random time-series')


%% correlation in dataset
ix = [5,7]; reps = 5;  % total number in correlated groups = #ixs + (#ixs * reps)
r_corr = r; 
for i=1:numel(ix)
    poss = randperm(N);
    r_corr(:,poss(1:reps)) = repmat(r(:,ix(i)),1,reps); % copy that column to reps others
end
%r_corrN = (r_corr-mean(r_corr(:)))./std(r_corr(:));
Cxy = corrcoef(r_corr);

[Cs,Cn,ixRetain,sigR,noiseR] = deconstructCxy(Cxy,T); 

figure
subplot(3,3,1), imagesc(Cxy); title('Cxy')
subplot(3,3,2), imagesc(Cs); title('Signal'); colorbar
subplot(3,3,3), imagesc(Cn); title('Noise')

subplot(3,3,7), imagesc(Cxy(ixRetain,ixRetain)); title('Cxy retained')
subplot(3,3,8), imagesc(sigR); title('Signal, retained'); colorbar
subplot(3,3,9), imagesc(noiseR); title('Noise, retained')

% get eigenvalues
[V,D] = eig(Cxy);
egs =diag(D);

% upper limit of random matrix distribution
MP_lambda_pos = (1+sqrt(N/T)).^2;
MP_lambda_neg = (1-sqrt(N/T)).^2;

subplot(3,3,[4:6])
hist(egs,50)
line([MP_lambda_neg MP_lambda_neg],[0 10],'Color',[1 0 0])
line([MP_lambda_pos MP_lambda_pos],[0 10],'Color',[1 0 0])
xlabel('Eigenvalues')
ylabel('Frequency')
title('Eigenvalue distribution of Cxy of random time-series with 2 correlated groups')


