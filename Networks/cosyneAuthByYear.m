%% Load cosyne authorship data and strip out author matrix for each year
% cosyneData is a structure that has:
% lastnames - a cell array of each occurance of extracted last names.
% firstnames - equivalent
% posterID - a cell array of posters with a semi-arbitrary ID. The format is IDNUMBER.YEAR, so 123.2011 means it is poster 123 from 2011
% authorID - an arbitrary ID number I assigned to each author. I forget how I came up with it.
% posterHash - not relevant here, but I also assign poster number to poster ID for an analysis that didn't pan out
% authorHash - for each author ID, the associated name.
    
clear all; close all;

load('cosyneFinalData');
A_Full = adjMatrix;

% get complete list of years
nP = cellfun('length',cosyneData.posterID);

years = fix(rem([cosyneData.posterID{:}],1)*10^4);  % just the year (after decimal)
IDs = floor([cosyneData.posterID{:}]);  % just the poster ID
yrSet = unique(years); 

% posterIDs can have multiple entries per author
% are these then posters from same author in same year? YES
% but different years treated as different authors? YES

% SO....
% for each year: find all poster cells first
% then expand into list of [authorID posterID]; find all matches as below


% loop over year list and create new matrices per year
for iY = 1:numel(yrSet)

    % get all poster IDs for that year
    ix = find(years == yrSet(iY));

    % get all corresponding author IDs for that year
    ixAuthor = cosyneData.authorID(ix);  % each author of each poster
    setA = unique(ixAuthor); % set of unique authors in that list
    
    % get all poster IDs
    ixPoster = IDs(ix);
    thesePosters = unique(ixPoster);
    
    % find out which were on the same posters...
    adjMatrix = zeros(numel(setA));
    for iP = 1:numel(thesePosters)
        matchPoster = find(ixPoster == thesePosters(iP));  % find all authors on this poster
        auths = ixAuthor(matchPoster);
        ixMatrix =arrayfun(@(x) find(x == setA),auths);
        adjMatrix(ixMatrix,ixMatrix) = adjMatrix(ixMatrix,ixMatrix) + 1;  %
    end
    
    % extract subset of names (firstname, lastname)
    nodelabels = cell(numel(setA),1);
    for iA = 1:numel(setA)
        iNames = ix(find(ixAuthor == setA(iA),1));  % get first appearance of this author in the current year, 
        % and get names from master list
        nodelabels{iA} = [cosyneData.firstnames{iNames}, cosyneData.lastnames{iNames}];
    end
    % save 
    fname = ['CosyneYear' num2str(yrSet(iY))];
    save(fname,'adjMatrix','nodelabels');
end