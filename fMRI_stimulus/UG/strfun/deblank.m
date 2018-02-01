%DEBLANK Remove trailing blanks
%   R = DEBLANK(S) removes any trailing whitespace and null characters from 
%   the end of S. However, DEBLANK does not remove significant whitespace
%   characters (such as (char(160)).
%
%   S can be a character vector, cell array of character vectors, or a
%   string array. If S is an array, then DEBLANK removes trailing blanks from 
%   each element. DEBLANK returns R as the same data type as S.
%
%   Examples:
%   A{1,1} = 'MATLAB    ';
%   A{1,2} = 'SIMULINK    ';
%   A = deblank(A)
%   A = 
%      'MATLAB'    'SIMULINK'
%       
%   See also ISSPACE, CELLSTR, STRING, STRTRIM.

%   Copyright 1984-2016 The MathWorks, Inc.
%==============================================================================

