function varargout = split(str, varargin)
%SPLIT Split strings in string array
%   NEWSTR = SPLIT(STR) splits STR at whitespace into string array NEWSTR.
%   If STR is a string scalar or character vector, then NEWSTR is a P-by-1
%   string array where P is the number of splits.  Whitespace is defined as
%   any sequence of whitespace characters such as spaces, tabs, and
%   newlines.
%
%   STR also can be a character vector or a cell array of character
%   vectors.
% 
%   If STR is an M-by-1 string vector, then the size of NEWSTR is M-by-P.
%   If STR is 1-by-N, then NEWSTR is 1-by-N-by-P. If STR is a strings
%   M-by-N-by-... array, NEWSTR is a string array with an additional
%   dimension containing the number of splits P from each element in STR.
%   The size of NEWSTR is M-by-N-by-...-by-P.
% 
%   If the number of splits P is not the same for every element in STR,
%   then SPLIT will error.  In that case, write a FOR loop to call SPLIT on
%   each element of STR.
% 
%   NEWSTR = SPLIT(STR,DELIMITER) splits STR at DELIMITER into NEWSTR.
%   DELIMITER can be a string array, a character vector, or a cell array of
%   character vectors. If DELIMITER is a string array or a cell array,
%   SPLIT splits STR along the elements in DELIMITER, in the order in which
%   they appear in the DELIMITER array.
%
%   NEWSTR = SPLIT(STR,DELIMITER,DIM) splits STR at DELIMITER into NEWSTR
%   along the specified dimension. The default value of DIM is the first
%   trailing dimension equal to 1.
% 
%   [NEWSTR,MATCHES] = SPLIT(...) returns the string array, MATCHES,
%   containing all occurrences of the delimiters at which STR was split.
%   MATCHES always contains one fewer element than NEWSTR along the
%   dimension P.
%
%   Example:
%   names = string({'Mary Jones';'John Smith';'Elizabeth Young'});
%   names = split(names)
%
%   returns
%
%       "Mary"         "Jones"
%       "John"         "Smith"
%       "Elizabeth"    "Young"
%
%   Example:
%   names = string({'Mary Jones';'John Smith';'Elizabeth Young'});
%   names = split(names,' ',1)
%
%   returns
%
%       "Mary"     "John"     "Elizabeth"
%       "Jones"    "Smith"    "Young"    
%
%   Example:
%   myPath = string('/Users/jdoe/My Documents/Examples');
%   myFolders = split(myPath,'/')
%
%   returns
%
%       ""
%       "Users"
%       "jdoe"
%       "My Documents"
%       "Examples"
%
%   See also SPLITLINES, JOIN, COMPOSE, NEWLINE, COUNT, STRING

%   Copyright 2015-2016 The MathWorks, Inc.

    narginchk(1, 3);
    if ~(ischar(str) && (isempty(str) || isrow(str))) && ~iscellstr(str) && ~isstring(str)
        error(fillMessageHoles('MATLAB:string:MustBeCharCellArrayOrString',...
                               'MATLAB:string:FirstInput'));
    end

    try
        stringStr = string(str);
        [varargout{1:nargout}] = stringStr.split(varargin{:});
    catch E
        throw(E)
    end
end
