function s = insertAfter(str, pos, text)
% INSERTAFTER Insert substring after a specified position
%    S = INSERTAFTER(STR, START_STR, NEW_TEXT) returns a string with the
%    content of STR where NEW_TEXT is inserted directly after every
%    occurrence of START_STR. S has the same dimensions as STR. START_STR
%    can be a string scalar or the same size as STR.
%
%    S = INSERTAFTER(STR, POS, NEW_TEXT) returns a string with the content
%    of STR where NEW_TEXT is inserted directly after numeric position POS.
%    POS must be an interger from 1 to the number of characters in a string
%    element.
%
%    NOTE: STR, START_STR, and NEW_TEXT can also be character vectors or
%    cell arrays of character vectors. The output S is the same data type
%    as STR.
%
%    Example:
%        str = string('The quick fox');
%        insertAfter(str,'quick',' brown')
%
%        returns 
%
%            "The quick brown fox"
%
%    See also STRING/PLUS, INSERTBEFORE, REPLACE, REPLACEBETWEEN

%   Copyright 2016 The MathWorks, Inc.

    narginchk(3, 3);
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
        s = convert(stringStr.insertAfter(pos, text));
    catch ex
        if strcmp(ex.identifier, 'MATLAB:string:CannotConvertMissingElementToChar')
            error(fillMessageHoles('MATLAB:string:CannotInsertMissingIntoChar',...
                                   'MATLAB:string:MissingDisplayText'));
        end
        ex.throw;
    end
end
