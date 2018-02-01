function s = extractBetween(str, start, stop, varargin)
% EXTRACTBETWEEN Extract substring from text
%    S = EXTRACTBETWEEN(STR, START_STR, END_STR) returns the substring of
%    STR that occurs between the strings START_STR and END_STR but does not
%    include either START_STR and END_STR. If multiple non-overlapping
%    START_STR END_STR pairs are found, S grows along the first trailing
%    dimension whose size is 1 and all the pairs are returned. START_STR
%    and END_STR can be string scalar the same size as STR.
% 
%    S = EXTRACTBETWEEN(STR, START, END) returns a substring of STR that
%    starts at numeric position START and ends at numeric position END.
%    START and END must be intergers between and including 1 to the number
%    of characters in a string element. The substring S includes the
%    characters at those positions.
% 
%    NOTES: 
%    * STR, START_STR, and END_STR can also be character vectors or
%      cell arrays of character vectors. The output S is always a string.
%
%    * The START_STR and END_STR strings are an open set (i.e., bounds
%      are not included in the output) while the START and END numeric
%      positions are a closed set (i.e., bounds are included in the output).
% 
%    * Position and string bounds can be mixed.  Their inclusivity persists
%      with the argument type.
%
%    S = EXTRACTBETWEEN(..., 'Boundaries', B) forces the bounds to be
%    inclusive when B is 'inclusive' and forces the bounds to be exclusive
%    when B is 'exclusive'.
%
%    Example:
%        str = string('The quick brown fox');
%        extractBetween(str,'quick ',' fox')
%
%        returns 
%
%            "brown"
%
%    See also ERASEBETWEEN, REPLACEBETWEEN, EXTRACTBEFORE, EXTRACTAFTER

%   Copyright 2016 The MathWorks, Inc.

    narginchk(3, Inf);
    if ~(ischar(str) && (isempty(str) || isrow(str))) && ~iscellstr(str) && ~isstring(str)
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    end

    try
        stringStr = string(str);
        s = stringStr.extractBetween(start, stop, varargin{:});
    catch E
        throw(E);
    end
end
