% script to unit test the conversion to multiple edges
% Mark Humphries 7/9/2020
clearvars

% parameters of synthetic network
N = 100;
density = 0.2;
binary_threshold = 0.1; % everything below this set to zero
integer_conversion = 10;  % convert to integer edge count

% cases to check
C_to_check = {-1,0.7,'alll','all',2,10,100};
 
%% make synthetic weight matrices
W = triu(rand(N),1);
idxs = find(triu(ones(N),1));  % all indexes of upper triangular
zeroed = idxs(randperm(numel(idxs)));       % pick (1-density) of them to be zero 
zeroed = zeroed(1:round((1-density)*numel(idxs)));  % number of upper triangular entries
W(zeroed) = 0;
W = W + W';               % make undirected, real-valued network

% find minimum weight value
minW = min(min(W(W>0)));

% make binary version
W_binary = W > binary_threshold;

% make integer version
W_int = round(W*integer_conversion);

% expected answers for each C checked
W_expected = {[],[],[],round(W*1/minW),round(W*C_to_check{5}),round(W*C_to_check{6}),round(W*C_to_check{7})};
W_binary_expected = {[],[],[],W_binary,W_binary,W_binary,W_binary};
W_int_expected = {[],[],[],W_int,W_int,W_int,W_int};

%% check weight matrix

nZero = sum(sum(W==0)) - N;             % count number of non-zero entries, correct count for diagonal 
Wdensity = 1 - (nZero / (N*(N-1)));     % density of weight matrix


%%  check the conversion function itself
for iC = 1:numel(C_to_check)    
    pars = {W,C_to_check{iC}};
    [blnExpectedFcn,A,conversion_used] = doUnitTest('ConvertToMultiEdges',pars);
    % check conversion returned by function is correct (as it's used later
    % in the main functions to reconvert the generated null model networks back to the
    % same value range as the data network)
    if ~blnExpectedFcn
        conversion_used
    end
    % check each A to make sure it has correct conversion
    if nnz(A ~= W_expected{iC}) > 0
        error(['Returned converted matrix for Real valued W is not correct for checked input ' num2str(iC) ])
    end
    
    % binary networks
    pars = {W_binary,C_to_check{iC}};
    [blnExpectedFcn(iC),A] = doUnitTest('ConvertToMultiEdges',pars);
    % check each A to make sure it has correct conversion
    if nnz(A ~= W_binary_expected{iC}) > 0
        error(['Returned converted matrix for Binary W is not correct for checked input ' num2str(iC) ])
    end
    
    % integer networks
    pars = {W_int,C_to_check{iC}};
    [blnExpectedFcn(iC),A] = doUnitTest('ConvertToMultiEdges',pars);
    % check each A to make sure it has correct conversion
    if nnz(A ~= W_int_expected{iC}) > 0
        error(['Returned converted matrix for Integer W is not correct for checked input ' num2str(iC) ])
    end

end


%% check calling it in the scripts
% not checking answers, just checking that each call is correctly
% implemented in code
for iC = 1:numel(C_to_check)  
    % real-valued networks
    try
       E = poissonSparseWCM(W,10,C_to_check{iC});
       E = poissonFullWCM(W,10,C_to_check{iC});
       E = linkSparseWCM(W,10,C_to_check{iC});
       E = linkFullWCM(W,10,C_to_check{iC});
       
    catch ME
        disp(ME.identifier) % custom error checking to deal with nested errors
        if strfind(ME.identifier,'MATLAB')
            rethrow(ME)
        end
    end
end

