function s = erase(str, match)
%ERASE Remove content from text
%   MODIFIEDSTR = ERASE(ORIGSTR,MATCHSTR) removes all
%   occurrences of MATCHSTR from ORIGSTR.
%
%   ORIGSTR and MATCHSTR can be string, char, or cell arrays of 
%   character vectors. MATCHSTR can be any sized array. Each occurrence 
%   of all the elements of MATCHSTR will be removed from ORIGSTR.
% 
%   See also STRREP, REGEXPREP, REPLACE, ERASEBETWEEN

%   Copyright 2016 The MathWorks, Inc.

    narginchk(2, 2);
    if ischar(str) && (isempty(str) || isrow(str))
        convert = @char;
    elseif iscellstr(str)
        convert = @cellstr;
    elseif isstring(str)
        convert = @string;
    else
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    end

    try
        stringStr = string(str);
        s = convert(stringStr.erase(match));
    catch E
        throw(E)
    end
end
