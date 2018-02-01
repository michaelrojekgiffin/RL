function s = join(str, varargin)
% JOIN Append elements of a string array together
%    NEWSTR = JOIN(STR) appends the elements of STR, placing a space
%    character between consecutive strings, and returns the result as a
%    string array NEWSTR. JOIN will combine the strings along the last
%    dimension of STR not equal to 1. STR can be a string array or cell
%    array of character vectors.
% 
%    NEWSTR = JOIN(STR,DELIMITER) appends the elements of STR and places
%    elements of DELIMITER between them. If DELIMITER is an array, then it
%    must have one element less than STR along the dimension being joined.
%    The size of every other dimension of DELIMITER either must be 1 or
%    must match the size of the corresponding dimension of STR. The space
%    character is the default value of DELIMITER.
% 
%    NEWSTR = JOIN(STR,DIM) appends the elements of STR along the dimension
%    DIM. The default value of DIM is the last dimension of STR with a size
%    that does not equal 1.
% 
%    NEWSTR = JOIN(STR,DELIMITER,DIM) appends the elements of STR along
%    the dimension DIM and places elements of DELIMITER between the
%    strings.
%
%    Example:
%        STR = string({'John','Smith';'Mary','Jones'});
%        join(STR)
%    
%    returns
%
%        "John Smith"
%        "Mary Jones"
%
%    Example:
%        STR = string({'John','Smith';'Mary','Jones'});
%        join(STR,1)
%
%    returns
%
%        "John Mary"    "Smith Jones"
%
%    Example:
%        STR = string({'x','y','z';'a','b','c'});
%        DELIMITER = {' + ',' = ';' - ',' = '};
%        join(STR,DELIMITER)
%
%    returns
%
%        "x + y = z"
%        "a - b = c"    
%
%    See also SPLIT, STRING/PLUS, COMPOSE

%   Copyright 2015-2016 The MathWorks, Inc.

    narginchk(1, 3);
    if ~iscellstr(str) && ~isstring(str)
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeStringOrCellArrayOfStrings', firstInput));
    end

    try
        stringStr = string(str);
        s = stringStr.join(varargin{:});
    catch E
        throw(E);
    end
end
