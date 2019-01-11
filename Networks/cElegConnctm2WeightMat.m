%% script to compute weight matrices of C. elegans neuron network
%
% 11/11/2018: created
%
% Javier Caballero, 

clc
clear



%% import data
filename = 'cElegansConnectome.csv';
delimiter = ',';
startRow = 2;

% format per line
formatSpec = '%C%C%C%f%[^\n\r]';

% open file
fileID = fopen(filename,'r');

% read columns of data according to the format
data = textscan(fileID, formatSpec, 'Delimiter', ...
    delimiter, 'TextType', 'string', 'EmptyValue', NaN, ...
    'HeaderLines' ,startRow-1, 'ReturnOnError', false, ...
    'EndOfLine', '\r\n');

% close file
fclose(fileID);

% make into table
data = table(data{1:end-1}, 'VariableNames', ...
    {'Neuron1','Neuron2','Type','Nbr'});



%% make weight matrices
% verify send + received chemical synapse integrity (S + Sp = R + Rp)
sentEqual2Received = (sum(data.Nbr(data.Type == 'S' | ...
    data.Type == 'Sp')) - sum(data.Nbr(data.Type == 'R' |...
    data.Type == 'Rp'))) == 0

% unique names of neurons
neuronNames = sort(unique(data.Neuron1));

% types of synapse...
% ... all
typesAll = unique(data.Type);
% remove neuro-muscular junctions and redundant synapses
typesAll(typesAll == 'NMJ' | typesAll == 'R' | typesAll == 'Rp') = [];
% ... all chemical
typesChem = categorical({'S'; 'Sp'});
% ... all electric
typesElect = 'EJ';

% fill weight matrices
for countMat = 1:3% matrix-type-wise
    
    % initialize weight matrix
    wMat = zeros(size(neuronNames, 1), size(neuronNames, 1));
    
    % types of synapse to be used
    if countMat == 1
        types = typesAll;
    elseif countMat == 2
        types = typesChem;
    else
        types = typesElect;
    end
    
    % connection-wise
    for countSend = 1:size(neuronNames, 1)% sending neuron
        dataTemp = data(ismember(data.Type, types) & ...
            data.Neuron1 == neuronNames(countSend), :);
        
        for countRec = 1:size(dataTemp, 1)% receiving neuron
            wMat(countSend, neuronNames == ...
                dataTemp.Neuron2(countRec)) = wMat(countSend,...
                neuronNames == dataTemp.Neuron2(countRec)) +...
                dataTemp.Nbr(countRec);
        end
        
    end
    
    % name martrix
    if countMat == 1
        wMatAll = wMat;
    elseif countMat == 2
        wMatChem = wMat;
    else
        wMatElect = wMat;
    end
    
    %  plot
    figure(countMat)
    imagesc(wMat)
    colorbar
    
    % set colormap
    colormap(flipud(bone))
end

% make undirected versions
wMatAllUndirected = wMatAll + wMatAll' - wMatElect;% electric synapses are taken as single links
wMatChemUndirected = wMatChem + wMatChem';
wMatElectUndirected = wMatElect;% this was already undirected



%% save
% weight matrices as directed networks
save('cElegWeightMatAllSynap.mat', 'wMatAll')% all synapses
save('cElegWeightMatChemSynap.mat', 'wMatChem')% only chemical
save('cElegWeightMatElectSynap.mat', 'wMatElect')% only electrical
% weight matrices as undirected networks
save('cElegWeightMatAllSynapUndirected.mat', 'wMatAllUndirected')% all synapses
save('cElegWeightMatChemSynapUndirected.mat', 'wMatChemUndirected')% only chemical
save('cElegWeightMatElectSynapUndirected.mat', 'wMatElectUndirected')% only electrical
% ordered names of neurons
save('cElegNeuronList.mat', 'neuronNames')


