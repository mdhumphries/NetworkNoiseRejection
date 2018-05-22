function [blnExpectedError,varargout] = doUnitTest(fcn,pars)

% DOUNITTEST runs the specified function and deals with errors
% [B,V] = DOUNITTEST(F,P) runs the function specified by the string F,
% given the arguments in the cell array P
%
% Returns:
%   B: boolean flag (0,1) indicating whether an expected error was thrown
%   or not
%   V: the set of output arguments from function F; V is a cell array of
%   the length of output arguments expected from function F (which is 
%   defined by the number of output arguments in the call to this function).
%   For example, calling [B,x1,x2] = feval(F,P) will expect 2 output
%   arguments from F, and return them in x1 and x2.
%       If B=1, then all arguments are empty ([]).
%
% IMPORTANT:
%   Within each function, the expected error messages must be identified by
%   the name of the function i.e. error(F,'some error text');
%   DOUNITTEST searches for the name of error message within the Message
%   Identifier to know that this was an expected error.
%
% Mark Humphries 7/3/2017

blnExpectedError = 0;
V = {};
nout = nargout-1;

try
    %V = feval(fcn,pars{:});
    [V{1:nout}] = feval(fcn,pars{:});
    % assign outputs according to whether the function executed or not

    for iN = 1:nout
        varargout(iN) = V(iN);
    end
catch ME
    if any(strfind(ME.identifier,fcn))
        strPars = [];
        for iP = 1:numel(pars)
            if ischar(pars{iP})
                strAdd = pars{iP};
            elseif iscell(pars{iP})
                strAdd = 'cell';
            elseif isstruct(pars{iP})
                strAdd = 'struct';
            elseif numel(pars{iP}) > 1
                d = size(pars{iP}); strDim = [];
                for iD = 1:numel(d)
                    strDim = [strDim num2str(d(iD)) '*'];
                end
                strAdd = [strDim(1:end-1) ' array'];
            else
                strAdd = num2str(pars{iP});
            end               
            strPars = [strPars 'Par(' num2str(iP) ') = ' strAdd '\n'];  
        end
        msg = sprintf(['\n ' fcn ' Error: \n' ME.message ' \n Thrown by: \n ' strPars]); 
        disp(msg)
        blnExpectedError = 1;
        varargout = cell(nout,1);
    else
       rethrow(ME); 
    end
end


