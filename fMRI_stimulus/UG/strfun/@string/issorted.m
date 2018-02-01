%ISSORTED Determine whether string elements are in sorted order
%   TF = ISSORTED(STR) returns logical 1 (true) if the elements of string
%   array STR are sorted in ascending order and logical 0 (false) otherwise.
%
%   MATLAB stores characters as Unicode using the UTF-16 character encoding scheme. 
%   SORT sorts strings according to the UTF-16 code point order. For the characters 
%   that are also the ASCII characters, this order means that uppercase letters come before 
%   lowercase letters. Digits and many punctuation marks also come before letters.
%
%   If STR is an M-by-N string array, ISSORTED returns true if each column is
%   sorted. If STR is a multidimensional array, ISSORTED returns true if STR is
%   sorted along the first array dimension whose size does not equal 1.
%
%   TF = ISSORTED(S,DIM) returns true if STR is sorted along dimension DIM.
%   For N dimensions, DIM is an integer between 1 and N.
%
%   TF = ISSORTED(S,DIM,DIRECTION) returns true if STR is sorted along the
%   direction specified by DIRECTION. DIRECTION is 'ascend' for ascending order, or 
%   'descend' for descending order. The default value of DIRECTION is
%   'ascend'.
%
%   Example:
%       STR = string({'Jones','Evans'});
%       issorted(STR)                     
%
%       returns 
%
%          0
%
%       issorted(S,'descend')             
%
%       returns 
%
%          1
%
%   Example:
%       STR = string({'Jones','Evans';'Smith','Peterson'});
%
%       issorted(STR)                     
%
%       returns  
%
%          1
%
%       issorted(STR,2)                   
%
%       returns  
%
%          0
%
%       issorted(STR,2,'descend')         
%
%       returns  
%
%          1
%
%   See also STRING/SORT.

%   Copyright 2015-2016 The MathWorks, Inc.

