function s = extractAfter(str, pos)
% EXTRACTAFTER Extract substring after a specified position
%    S = EXTRACTAFTER(STR, START_STR) returns a substring of STR that
%    starts with the character after the first occurrence of START_STR and
%    ends with the last character of STR. If STR is an array, S will be an
%    array with the same dimensions as STR containing the substrings.
%    START_STR can be a string scalar or an array the same size as STR.
%
%    S = EXTRACTAFTER(STR, POS) returns a substring of STR that that starts
%    with the character after numeric position POS and ends with the last
%    character of STR. POS must be an interger from 1 to the number of
%    characters in a string element.
%
%    NOTE: STR and START_STR can also be character vectors or cell arrays
%    of character vectors. The output S is always a string.
%
%    Example:
%        str = string('The quick brown fox');
%        extractAfter(str,'quick ')
%
%        returns 
%
%            "brown fox"
%
%    See also EXTRACTBEFORE, EXTRACTBETWEEN

%   Copyright 2016 The MathWorks, Inc.

    narginchk(2, 2);
    if ~(ischar(str) && (isempty(str) || isrow(str))) && ~iscellstr(str) && ~isstring(str)
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    end

    try
        stringStr = string(str);
        s = stringStr.extractAfter(pos);
    catch E
        throw(E);
    end
end
