function s = extractBefore(str, pos)
% EXTRACTBEFORE Extract substring before a specified position
%    S = EXTRACTBEFORE(STR, END_STR) returns a substring of STR that starts
%    with the first character of STR and ends with the last character
%    before the first occurrence of END_STR. If STR is an array, S will be
%    an array with the same dimensions as STR containing the substrings.
%    END_STR can be a string scalar or an array the same size as STR.
%
%    S = EXTRACTBEFORE(STR, POS) returns a substring of STR that begins
%    with the first character of STR and ends with the last character
%    before the numeric position specified by POS. POS must be an interger
%    from 1 to the number of characters in a string element.
%
%    NOTE: STR and END_STR can also be character vectors or cell arrays of
%    character vectors. The output S is always a string.
%
%    Example:
%        str = string('The quick brown fox');
%        extractBefore(str,' brown')
%
%        returns 
%
%            "The quick"
%
%    See also EXTRACTAFTER, EXTRACTBETWEEN

%   Copyright 2016 The MathWorks, Inc.

    narginchk(2, 2);
    if ~(ischar(str) && (isempty(str) || isrow(str))) && ~iscellstr(str) && ~isstring(str)
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    end

    try
        stringStr = string(str);
        s = stringStr.extractBefore(pos);
    catch E
        throw(E);
    end
end
