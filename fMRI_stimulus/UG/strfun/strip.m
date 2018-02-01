function s = strip(m, varargin)
%STRIP Remove leading and trailing whitespaces
%   NEWSTR = STRIP(STR) removes all consecutive whitespace characters from
%   the beginning and the end of each element in STR. Whitespace is defined
%   as any sequence of whitespace characters such as spaces, tabs, and
%   newlines. 
%
%   NEWSTR = STRIP(STR,SIDE) removes whitespace characters from the
%   specified SIDE. SIDE can be 'left', 'right', or 'both'.  The default
%   value of SIDE is 'both'.
%
%   NEWSTR = STRIP(STR,PAD_CHARACTER) removes PAD_CHARACTER from STR.
%   PAD_CHARACTER must be exactly one character.
% 
%   NEWSTR = STRIP(STR,SIDE,PAD_CHARACTER) removes PAD_CHARACTER from the
%   specified SIDE.
% 
%   Example:
%
%       STR = string({'moustache '; 
%                     '   goatee';
%                     '   beard    '});
%       strip(STR)
%
%       returns
%
%           "moustache"
%           "goatee"
%           "beard"
%
%   Example:
%       strip(string('C:\Temp\Files\'),'right','\')
%
%       returns
%
%           "C:\Temp\Files"
%       
%   See also PAD, STRING, REPLACE

%   Copyright 2016 The MathWorks, Inc.

    narginchk(1, 3);
    if ischar(m) && (isempty(m) || isrow(m))
        convert = @char;
    elseif iscellstr(m)
        convert = @cellstr;
    elseif isstring(m)
        convert = @string;
    else
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    end

    try
        stringM = string(m);
        s = convert(stringM.strip(varargin{:}));
    catch E
        throw(E)
    end
end
