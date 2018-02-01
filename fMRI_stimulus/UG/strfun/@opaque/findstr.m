function s = findstr(s1,s2)
%FINDSTR Find one string within another for Java objects.
%
%   FINDSTR will be removed in a future release. Use STRFIND instead.

%   Copyright 1984-2009 The MathWorks, Inc.

s = findstr(fromOpaque(s1),fromOpaque(s2));




