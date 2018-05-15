%% testing Marchenko-Pastur for actual random matrices


% vrX = var(X(:));
% Y = 1/N * X * X';  % normalised matrix
% [V,D] = eig(Y,'vector');
% 
% maxEg = vrX * (1+sqrt(c)).^2;  % upper M-P bound
% 
% histogram(D);
% line([maxEg maxEg],[0 10])

%% Wigner
% Gaussian
N = 1000; c=1;
X = randn(N);  % asymmetric matrix

H = (X + X') / 2;
Es = eig(H);
En = Es / sqrt(N/2);
dx = 0.05;
x = -2:dx:2;

counts = hist(En,x); 
W = 1/(2*pi) * sqrt(4-x.^2);

figure
bar(x,counts / (dx * N)); hold on
plot(x,W,'r')


%% Tracey-Widom for complex numbers
Emax = []; rpts = 1000; N = 100;
for i=1:rpts
    X = randn(N) +sqrt(-1)*randn(N);  % complex numbers
    H = (X + X') / 2;
    Es = eig(H);
    Emax = [Emax; max(Es)];
end
dx = 0.01;
x = -5:dx:5;
EMaxScale = N.^(1/6)*(Emax - 2*sqrt(N)); % scale values; for f2 (complex numbers)

counts = hist(EMaxScale,x); 

[t,f2] = TraceyWidom(-5,5,dx);

figure
bar(x,counts/(rpts*dx)); hold on
plot(t,f2,'r')


