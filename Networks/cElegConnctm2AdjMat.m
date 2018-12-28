clc
clear



%% import data
filename = 'D:\Kronos\Chamba\Postdoc\Manchester\Documentation\Publications\netNoiseRejection\cElegansNetwork\cElegansConnectome.csv';
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



%% make adjacency matrices
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

% fill adjacency matrices
for countMat = 1:3% matrix-type-wise
    
    % initialize adjacency matrix
    adjMat = zeros(size(neuronNames, 1), size(neuronNames, 1));
    
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
            adjMat(countSend, neuronNames == ...
                dataTemp.Neuron2(countRec)) = adjMat(countSend,...
                neuronNames == dataTemp.Neuron2(countRec)) +...
                dataTemp.Nbr(countRec);
        end
        
    end
    
    % make symmetric version
    % symAdjMAt = triangle(adjMat)
    
    % name martrix
    if countMat == 1
        adjMatAll = adjMat;
    elseif countMat == 2
        adjMatChem = adjMat;
    else
        adjMatElect = adjMat;
    end
    
    % set colormap
    colormap(flipud(bone))
    
    %  plot
    figure(countMat)
    imagesc(adjMat)
    colorbar
end



%% save
save('cElegAdjMatAllSynap.mat', 'adjMatAll')
save('cElegAdjMatChemSynap.mat', 'adjMatChem')
save('cElegAdjMatElectSynap.mat', 'adjMatElect')
save('cElegNeuronList.mat', 'neuronNames')


