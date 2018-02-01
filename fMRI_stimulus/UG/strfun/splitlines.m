function s = splitlines(str)
%SPLITLINES Split string at newline characters
%   NEWSTR = SPLITLINES(STR) splits the string STR at newline characters into 
%   a string array NEWSTR.
% 
%   NEWSTR = SPLITLINES(STR) is equivalent to calling SPLIT and specifying
%   newline characters as the delimiters to split upon, as in the code:
%       NEWSTR = split(STR,compose(string({'\r\n', '\n', '\r'})))
%
%   If STR contains literals such as '\n' to indicate newlines, then
%   convert the literals to actual newline characters with the COMPOSE
%   or SPRINTF functions. You also can add newline characters to strings with the
%   NEWLINE function.
%
%   Example:
%       STR = string('Name:');
%       STR = STR + '\n' + 'Title:\n' + 'Company:';
%       STR = compose(STR);
%       splitlines(STR)
%
%       returns
%
%           "Name:"
%           "Title:"
%           "Company:"
%
%   Example:
%       STR = string('In Xanadu did Kubla Khan');
%       STR = STR + newline + 'A stately pleasure-dome decree';
%       splitlines(STR)
%
%       returns
%
%           "In Xanadu did Kubla Khan"
%           "A stately pleasure-dome decree"
%
%   See also SPLIT, COMPOSE, NEWLINE, SPRINTF.

%  Copyright 2015-2016 The MathWorks, Inc.

    narginchk(1, 1);
    if ~(ischar(str) && (isempty(str) || isrow(str))) && ~iscellstr(str)
        error(fillMessageHoles('MATLAB:string:MustBeCharCellArrayOrString',...
                               'MATLAB:string:FirstInput'));
    end

    try
        stringStr = string(str);
        s = stringStr.splitlines;
    catch E
       throw(E); 
    end
end
