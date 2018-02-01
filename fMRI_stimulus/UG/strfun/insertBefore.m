function s = insertBefore(str, pos, text)
% INSERTBEFORE Insert substring before a specified position
%    S = INSERTBEFORE(STR, END_STR, NEW_TEXT) returns a string with the
%    content of STR where NEW_TEXT is inserted directly before every
%    occurrence of END_STR. S has the same dimensions as STR. END_STR can
%    be a string scalar or the same size as STR.
%
%    S = INSERTBEFORE(STR, POS, NEW_TEXT) returns a string with the content
%    of STR where NEW_TEXT is inserted directly before numeric position
%    POS. POS must be an interger from 1 to the number of characters in a
%    string element.
%
%    NOTE: STR, END_STR, and NEW_TEXT can also be character vectors or cell
%    arrays of character vectors. The output S is the same data type as
%    STR.
%
%    Example:
%        str = string('The quick fox jumps');
%        insertBefore(str,' fox',' brown')
%
%        returns 
%
%            "The quick brown fox jumps"
%
%    See also INSERTAFTER, REPLACE, REPLACEBETWEEN, EXTRACTBEFORE,
%    EXTRACTAFTER, EXTRACTBETWEEN

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
        s = convert(stringStr.insertBefore(pos, text));
    catch ex
        if strcmp(ex.identifier, 'MATLAB:string:CannotConvertMissingElementToChar')
            error(fillMessageHoles('MATLAB:string:CannotInsertMissingIntoChar',...
                                   'MATLAB:string:MissingDisplayText'));
        end
        ex.throw;
    end

end
