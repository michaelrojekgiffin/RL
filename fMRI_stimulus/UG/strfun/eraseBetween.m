function s = eraseBetween(str, start, stop, varargin)
% ERASEBETWEEN Remove substring from text
%    S = ERASEBETWEEN(STR, START_STR, END_STR) removes the substring of STR
%    that occurs between the strings START_STR and END_STR but does not
%    remove START_STR and END_STR themselves. If multiple non-overlapping
%    START_STR and END_STR pairs are found, the characters between each
%    pair will be removed. START_STR and END_STR can be string scalar or
%    the same size as STR.
% 
%    S = ERASEBETWEEN(STR, START, END) removes the substring of STR from
%    numeric position START to numeric position END, including the
%    characters at those positions. START and END must be intergers between
%    and including 1 to the number of characters in a string element.
% 
%    NOTES: 
%    * STR, START_STR, and END_STR can also be character vectors or cell
%      arrays of character vectors. The output S is the same data type as
%      STR.
%
%    * The START_STR and END_STR strings are an open set (i.e., bounds are
%      included in the output) while the START and END numeric positions
%      are a closed set (i.e., bounds are not included in the output).
% 
%    * Position and string bounds can be mixed.  Their inclusivity persists
%      with the argument type.
%
%    S = ERASEBETWEEN(..., 'Boundaries', B) forces the bounds to be
%    inclusive when B is 'inclusive' and forces the bounds to be exclusive
%    when B is 'exclusive'.
%
%    Example:
%        str = string('The quick brown fox');
%        eraseBetween(str,'quick',' fox')
%
%        returns 
%
%            "The quick fox"
%
%    See also ERASE, EXTRACTBETWEEN, REPLACEBETWEEN, EXTRACTBEFORE,
%    EXTRACTAFTER

%   Copyright 2016 The MathWorks, Inc.

    narginchk(3, Inf);
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
        s = convert(stringStr.eraseBetween(start, stop, varargin{:}));
    catch E
        throw(E);
    end
end
