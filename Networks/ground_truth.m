%% Zacharay karate data set
load('karate.mat')
groups([1,2,3,4,5,6,7,8,11,12,13,14,17,18,20,22])=0;
groups([9,10,15,16,19,21,23,24,25,26,27,28,29,30,31,32,33,34])=1;

%% Dolphins dataset
load('dolphins.mat')
groups=zeros(62,1);
groups([2     6     7     8    10    14    18    20    23    26    27    28    32    33    42    49    55    57    58    61])=1;

%% adjnouns dataset
load('adjnoun.mat')
groups=Problem.aux.nodevalue;

%% polblogs dataset
load('polblogs.mat')
groups=Problem.aux.nodevalue;