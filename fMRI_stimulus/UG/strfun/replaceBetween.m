function s = replaceBetween(str, start, stop, value, varargin)
% REPLACEBETWEEN Replace bounded substring in string elements
%    S = REPLACEBETWEEN(STR, START_STR, END_STR, NEW_TEXT) replaces the
%    substring in STR that occurs between START_STR and END_STR with
%    NEW_TEXT but does not replace START_STR and END_STR themselves. If
%    multiple non-overlapping START_STR and END_STR pairs are found, the
%    characters between each pair will be replaced. START_STR and END_STR
%    can be string scalars or the same size as STR.
% 
%    S = REPLACEBETWEEN(STR, START, END, NEW_TEXT) replaces the substrings
%    in STR that occur between numeric positions START and END, including
%    the characters at those positions, with NEW_TEXT. START and END must
%    be intergers between and including 1 to the number of characters in a
%    string element.
%
%    NOTES: 
%    * STR, START_STR, END_STR, and NEW_TEXT can also be character
%      vectors or cell arrays of character vectors. The output S is the
%      same data type as STR.
% 
%    * The START_STR and END_STR strings are an open set (i.e., bounds
%      are included in the output) while the START and END numeric positions
%      are a closed set (i.e., bounds are not included in the output).
%
%    * Position and string bounds can be mixed.  Their inclusivity persists
%      with the argument type.
%
%    S = REPLACEBETWEEN(..., 'Boundaries', B) forces the bounds to be
%    inclusive when B is 'Inclusive' and forces the bounds to be exclusive
%    when B is 'Exclusive'.
%
%
%    Example:
%        str = string('The quick brown fox jumps');
%        replaceBetween(str,'quick','fox', ' red ')
%
%        returns 
%
%            "The quick red fox jumps"
%
%    See also REPLACE, ERASE, ERASEBETWEEN, EXTRACTBETWEEN

%   Copyright 2016 The MathWorks, Inc.

    narginchk(4, Inf);
    if ischar(str) && (isempty(str) || isrow(str))
        convert = @char;
    elseif iscellstr(str)
        convert = @cellstr;
    elseif isstring(str)
        convert = @string;
    else
        error(fillMessageHoles('MATLAB:string:MustBeCharCellArrayOrString',...
                               'MATLAB:string:FirstInput'));
    end

    try
        stringStr = string(str);
        s = convert(stringStr.replaceBetween(start, stop, value, varargin{:}));
    catch ex
        if strcmp(ex.identifier, 'MATLAB:string:CannotConvertMissingElementToChar')
            error(fillMessageHoles('MATLAB:string:CannotInsertMissingIntoChar',...
                                   'MATLAB:string:MissingDisplayText'));
        end
        ex.throw;
    end
end
