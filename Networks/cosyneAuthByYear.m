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

% how many posters per author entry
nP = cellfun('length',cosyneData.posterID);

% get list of years
firstID = cellfun(@(x) x(1),cosyneData.posterID);  % the first poster ID entry for every author entry
years = cast(rem(firstID,1)*10^4,'int32');  % just the year (after decimal)
yrSet = unique(years); 

% posterIDs can have multiple entries per author
% are these then posters from same author in same year? YES
% but different years treated as different authors? YES

% SO....
% for each year: find all poster cells first
% then expand into list of [authorID posterID]; find all matches as below


% loop over year list and create new matrices per year
for iY = 1:numel(yrSet)

    % get all author IDs for that year
    ix = find(years == yrSet(iY));
    ixAuthor = cosyneData.authorID(ix);  % each author listed for that year: some may be listed twice
    
    setA = unique(ixAuthor); % set of unique authors in that list
    n = arrayfun(@(x) sum(x == ixAuthor),setA);  % check where some authors appear more than twice
    
     % get a list of unique poster-author IDs and 
    ixPoster = floor([cosyneData.posterID{ix}]);
    thesePosters = unique(ixPoster);
    
    % matching list of author IDs
    thisNP = nP(ix);
    ixPA = arrayfun(@(x,y) zeros(x,1)+y,thisNP,ixAuthor,'UniformOutput',false);
    ixPosterAuthor = vertcat(ixPA{:});
    
    % find out which were on the same posters...
    adjMatrix = zeros(numel(setA));
    for iP = 1:numel(thesePosters)
        matchPoster = find(ixPoster == thesePosters(iP));  % find all authors on this poster
        auths = ixPosterAuthor(matchPoster);
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
    
    % remove diagonal! Store separately as per-author posters...
    nPosters = diag(adjMatrix);
    adjMatrix(eye(numel(setA))==1) = 0;
    
    % save 
    fname = ['CosyneYear' num2str(yrSet(iY))];
    save(fname,'adjMatrix','nodelabels','nPosters');
end