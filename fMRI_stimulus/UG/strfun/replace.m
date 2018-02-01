function s = replace(str, old, new) 
%REPLACE Replace segments from string elements
%   NEWSTR = REPLACE(STR,OLD,NEW) replaces all occurrences of the string OLD 
%   in string array STR with the string NEW. 
%
%   STR must be a string. OLD and NEW can be strings, character vectors, or 
%   cell arrays of character vectors.  If OLD and NEW both contain more than 
%   one string, then they must have the same number of elements. All
%   nonoverlapping occurrences of each element of OLD in STR are replaced
%   by the corresponding element of NEW.
%
%   STR and OLD both contain a string with no characters (""), REPLACE does
%   not replace "" with the contents of NEW. 
%
%   Example:
%       claim = string('This is a good example');
%       new_claim = replace(claim,'good','great')
%
%       returns
%
%       This is a great example.
%
%   Example:
%       c_files = string({'c:\cookies.m'; ...
%                         'c:\candy.m';   ...
%                         'c:\calories.m'});
%       d_files = replace(c_files,'c:','d:')
%
%       returns
%
%           "d:\cookies.m"
%           "d:\candy.m"
%           "d:\calories.m"
%
%   Example: 
%       STR = string({'Submission Date: 11/29/15\r';
%                     'Acceptance Date: 1/20/16\r';
%                     'Contact: john.smith@example.com\r\n'});
%       OLD = {'\r\n','\r'};
%       NEW = '\n';
%       replace(STR,OLD,NEW)
%
%       returns
%
%           "Submission Date: 11/29/15\n"
%           "Acceptance Date: 1/20/16\n"
%           "Contact: john.smith@example.com\n"
%    
%   See also STRREP, REGEXP, REGEXPREP.

%   Copyright 2015-2016 The MathWorks, Inc. 

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
        s = convert(stringStr.replace(old, new));
    catch ex
        if strcmp(ex.identifier, 'MATLAB:string:CannotConvertMissingElementToChar')
            error(fillMessageHoles('MATLAB:string:CannotInsertMissingIntoChar',...
                                   'MATLAB:string:MissingDisplayText'));
        end
        ex.throw;
    end

end 
