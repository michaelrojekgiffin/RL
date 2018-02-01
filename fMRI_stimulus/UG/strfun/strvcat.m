function s=strvcat(varargin)
%STRVCAT Vertically concatenate character vectors.
%   S = STRVCAT(T1,T2,T3,..) forms the character matrix S containing the
%   text T1,T2,T3,... as rows.  Automatically pads each character vector
%   with spaces in order to form a valid matrix. Each text parameter, Ti,
%   can itself be a character matrix.  This allows the creation of
%   arbitrarily large character arrays.  Empty character arrays in the
%   input are ignored.
%
%   S = STRVCAT(C), when C is a cell array of character vectors, passes
%   each element of C as an input to STRVCAT.  Empty character vectors in
%   the input are ignored.
%
%   STRVCAT('Hello','Yes') is the same as ['Hello';'Yes  '] except
%   that the padding is done automatically.
%
%   STRVCAT is not recommended. Use CHAR instead.
%
%   See also STRING, PAD, CHAR, STRCAT

%   Copyright 1984-2016 The MathWorks, Inc.

numinput = nargin;
if numinput == 1 && iscellstr(varargin{1})
  varargin = (varargin{1});
end
% find the empty cells 
notempty = ~cellfun('isempty',varargin);
% vertically concatenate the non-empty cells.
s = char(varargin{notempty});


