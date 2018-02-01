%ISSTRPROP Check if string elements are of a specified category.
%   B = ISSTRPROP(TEXT, C) returns in B the logical array of the same shape as
%   TEXT, confirming whether the elements of TEXT are of category C. The type
%   of TEXT may be any of a cell array, a string array, a character array, or
%   any of the MATLAB numeric types. If TEXT is a cell array or nonscalar string
%   array, then B is a cell array of the same shape as TEXT. The classification
%   of the elements in TEXT is done according to the Unicode definition of the
%   specified category. That is, if the numeric value of an element in the input
%   array falls within the range that defines a Unicode character category, then
%   this element is classified as of that category. This classification of
%   characters is subject to the locale setting of MATLAB. 
%
%   B = ISSTRPROP(TEXT, C, 'ForceCellOutput', CELLOUTPUT) forces B to be a 
%   cell array when CELLOUTPUT is true.
%
%   INPUT ARGUMENTS
%
%   TEXT can be any of type char, int8, uint8, int16, uint16, int32, uint32,
%   int64, uint64, double, cell array, and string arrays. Except for strings,
%   a cell array may contain arrays of the aforementioned types.
%   
%   Numbers of type double are converted to int32 according to MATLAB rules of
%   double-to-int conversion. Numbers of type int64 and uint64 bigger than
%   int32(inf) saturate to int32(inf). 
%
%   Argument C must be from the following set:
%   'alpha'     : classify TEXT as in the alphabetic letter range
%   'alphanum'  : classify TEXT as in the alphanumeric range
%   'cntrl'     : classify TEXT as in the range of control characters, char(0:20).
%   'digit'     : classify TEXT as in the range of numeric digits
%   'graphic'   : classify TEXT as in the range of graphic characters. These are
%               all values that represent characters NOT of the set
%               {unassigned, space, line separator, paragraph separator, control
%               characters, Unicode format control characters, private
%               user-defined characters, Unicode surrogate characters, Unicode
%               other characters}.
%   'lower'     : classify TEXT as in the range of lowercase letters
%   'print'     : classify TEXT as in the range of graphic characters, plus
%               char(32).
%   'punct'     : classify TEXT as in the range of punctuation characters
%   'upper'     : classify TEXT as in the range of uppercase letters
%   'wspace'    : classify TEXT as in the range of whitespace characters; this
%               range includes the ANSI C definition of whitespace, 
%               {' ','\t','\n','\r','\v','\f'}, in addition to a number of
%               other Unicode characters. 
%   'xdigit'    : classify TEXT as in the range of valid hexadecimal digits
%
%   EXAMPLES
%
%   B = isstrprop('abc123efg','alpha') returns B  => [1 1 1 0 0 0 1 1 1]
%   B = isstrprop('abc123efg','digit') returns B  => [0 0 0 1 1 1 0 0 0]
%   B = isstrprop('abc123efg','xdigit') returns B => [1 1 1 1 1 1 1 1 0]
%   B = isstrprop([97 98 99 49 50 51 101 102 103],'digit') returns 
%       B => [0 0 0 1 1 1 0 0 0]
%   B = isstrprop(int8([97 98 99 49 50 51 101 102 103]),'digit') returns 
%       B => [0 0 0 1 1 1 0 0 0]
%   B = isstrprop(['abc123efg';'abc123efg'],'digit') returns
%       B => [0 0 0 1 1 1 0 0 0; 0 0 0 1 1 1 0 0 0]
%   B = isstrprop({'abc123efg';'abc123efg'},'digit') returns
%       B => {[0 0 0 1 1 1 0 0 0]; [0 0 0 1 1 1 0 0 0]}
%   B = isstrprop(sprintf('abc\n'),'wspace') returns B  => [0 0 0 1]
%
%   See also ISCHAR, ISSPACE.

%   Copyright 2003-2016 The MathWorks, Inc.
