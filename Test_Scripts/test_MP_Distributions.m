%% testing Marchenko-Pastur for actual random matrices

% Gaussian
N = 1000; c=1;
X = randn(N);  % asymmetric matrix
vrX = var(X(:));

Y = 1/N * X * X';  % normalised matrix
[V,D] = eig(Y,'vector');

maxEg = vrX * (1+sqrt(c)).^2;  % upper M-P bound

histogram(D);
line([maxEg maxEg],[0 10])

% Wigner
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

% Tracey-Widom 

