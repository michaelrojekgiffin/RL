function out = strmatch(str,strs,flag)
%STRMATCH Cell array based character vector matching.
%   Implementation of STRMATCH for cell arrays of character vectors.
%
%   STRMATCH will be removed in a future release. Use startsWith or
%   validatestring instead.
%
%   See also STRMATCH, startsWith, VALIDATESTRING

%   Copyright 1984-2016 The MathWorks, Inc.

narginchk(2,3);

if isempty(strs), out = []; return; end
if iscellstr(str), str = char(str); end
if iscellstr(strs), strs = char(strs); end

if ~ischar(str) || ~ischar(strs)
  error(message('MATLAB:strmatch:InvalidInput'));
end
if (nargin==3) && ~ischar(flag)
  error(message('MATLAB:strmatch:InvalidFlagInput'));
end

if nargin==2
  out = strmatch(str,strs);
else
  out = strmatch(str,strs,flag);
end

