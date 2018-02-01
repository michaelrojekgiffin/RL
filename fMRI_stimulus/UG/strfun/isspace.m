%ISSPACE True for whitespace characters
%   For a character array CHR, ISSPACE(CHR) returns an array the same size 
%   as CHR containing logical 1 (TRUE) where the elements of CHR are
%   Unicode-represented whitespace characters and logical 0 (FALSE) where
%   they are not.  
%
%   Whitespace characters for which ISSPACE returns TRUE include tab, line
%   feed, vertical tab, form feed, carriage return, and space, in addition
%   to a number of other Unicode characters. 
%
%   Example
%      isspace('  Find spa ces ')
%      Columns 1 through 13 
%         1   1   0   0   0   0   1   0   0   0   1   0   0
%      Columns 14 through 15 
%         0   1
%     
%   See also ISLETTER, ISSTRPROP.
 
%   Copyright 1984-2016 The MathWorks, Inc.
%   Built-in function.

