function varargout = ismember(A, B, varargin)
%ISMEMBER True for set member
%   Implementation of ISMEMBER for string arrays.

%   Copyright 2016 The MathWorks, Inc.

    narginchk(2, inf);

    if ~ischar(A) && ~iscellstr(A) && ~isstring(A)
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    elseif ~ischar(B) && ~iscellstr(B) && ~isstring(B)
        secondInput = getString(message('MATLAB:string:SecondInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', secondInput));
    end
    
    stringA = string(A);
    stringB = string(B);
    [varargout{1:nargout}] = matlab.internal.language.callOrdinaryFunction(mfilename, stringA, stringB, varargin{:});

end
