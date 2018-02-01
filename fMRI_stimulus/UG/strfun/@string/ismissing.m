%ISMISSING Determine whether string is missing
%   TF = ISMISSING(STR) returns a logical array the same size as the
%   string array STR, with logical 1 (true) for each element of STR that is a 
%   missing element, and logical 0 (false) otherwise.
%
%   Example:
%       STR = string(NaN);
%       ismissing(STR)  
%
%       returns  
%
%       1
%
%   Example:
%       STR = string('John');
%       STR(4) = 'Mary'
%
%       returns
%
%           "John"    <missing>    <missing>    "Mary"
%
%       ismissing(STR)
%
%       returns
%
%          0   1   1   0
%
%   See also EQ, NE, ISEMPTY, STRLENGTH.

%   Copyright 2016 The MathWorks, Inc.
