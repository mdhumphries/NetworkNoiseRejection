% to test multiwaySpectCommDet.m

clc
clear all

% % from AS hackathon
% % synthetic network with 4 communities
% load('javier_test_network.mat')
% % adjacency matrix
% A_adjMat = network.adjacency;
% A_adjMat = A_adjMat > 99.5;
% % ground truth group membership list
% g_iGroupsVerts = network.membership + 1;

% network science data tested in Zhang & Newman, 2015
load('netScience.mat')
% get only the largest connected component members
largestComp = [31,32,33,34,35,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,70,71,72,73,78,79,80,81,91,95,96,97,98,99,100,101,114,115,121,122,127,128,129,132,133,134,135,136,137,150,151,152,153,158,159,163,164,185,186,187,191,192,193,194,195,204,217,218,219,220,221,222,223,224,225,226,234,235,236,252,253,263,264,265,266,267,268,269,270,277,278,279,280,281,282,283,284,285,286,302,303,304,305,306,307,308,309,310,311,317,318,327,328,329,330,331,332,343,344,345,346,347,348,371,372,376,377,378,379,402,403,404,405,406,417,418,429,430,443,444,451,461,462,463,464,465,466,467,473,474,475,486,489,490,491,492,493,494,500,501,502,503,504,508,509,510,515,516,517,518,547,548,549,550,551,562,563,575,576,577,586,587,588,590,591,592,593,596,597,610,611,612,613,639,640,641,642,647,653,654,655,656,657,658,675,676,677,686,693,694,695,696,697,698,699,701,702,703,709,710,711,716,717,718,719,730,737,738,739,740,757,758,759,760,761,762,763,764,765,766,771,772,773,774,775,776,777,789,840,841,854,864,865,866,867,868,893,894,895,909,935,945,946,956,957,964,965,978,984,985,986,987,1006,1009,1022,1023,1024,1025,1026,1027,1028,1031,1040,1041,1042,1082,1083,1084,1087,1088,1089,1090,1092,1093,1122,1123,1124,1131,1136,1137,1138,1139,1146,1163,1164,1173,1178,1179,1181,1182,1183,1190,1191,1192,1196,1197,1198,1215,1216,1217,1218,1222,1229,1230,1240,1256,1264,1283,1284,1296,1313,1314,1315,1316,1317,1342,1343,1344,1345,1348,1349,1362,1363,1364,1385,1386,1390,1395,1396,1397,1398,1405,1406,1407,1408,1409,1414,1415,1416,1417,1418,1452,1453,1461,1462,1469,1470,1471,1482,1483,1530,1550,1551,1552,1557,1558,1559,1561,1562];
% adjacency matrix
A_adjMat = A_adjMat(:, largestComp);
A_adjMat = A_adjMat(largestComp, :);

% number of trials
noTrials = 30;

nGroups = 20;

for count = 1:noTrials
    if exist('g_iGroupsVerts', 'var') == 1
        [bestPartition] =...
            multiwaySpectCommDet(A_adjMat, nGroups, g_iGroupsVerts); 
    else
        [bestPartition,maxQPartition] =...
            multiwaySpectCommDet(A_adjMat,nGroups);
    end
     bestNoOfGroups(count) = max(bestPartition);
     QmaxNoOfGroups(count) = max(maxQPartition);
end

figure(1)
subplot(211),
hist(bestNoOfGroups, 1:max(bestNoOfGroups))
ylabel('frequency')
xlabel('No. of groups detected: knee')
subplot(212),
hist(QmaxNoOfGroups, 1:max(QmaxNoOfGroups))
ylabel('frequency')
xlabel('No. of groups detected: Qmax')

