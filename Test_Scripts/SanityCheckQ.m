
clear all; close all;

load UnitTestClusteringData
n = size(Data{1},1); 

perfectC{1} = randi(4,n,1);
perfectC{2} = [ones(n/2,1); ones(n/2,1)+1];
perfectC{3} = [ones(n/4,1); ones(n/4,1)+1; ones(n/4,1)+2; ones(n/4,1)+3];

for iD = 1:numel(Data)
    % get modularity: do this with configuration model to get appropriate data
    P{iD} = expectedA(Data{iD});
    B{iD} = Data{iD} - P{iD};  % also catch expectedA errors here
    m = sum(sum(Data{iD}))/2;  % sum of unique weights
    Q(iD) = computeQ(perfectC{iD},B{iD},m);
end

testC = grps{2};
m = sum(sum(Data{iD}))/2;
computeQ(testC,B{2},m);
