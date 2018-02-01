function t = strcat(varargin)
%STRCAT String horizontal concatenation.
%   T = STRCAT(S1,S2,...), when any of the inputs is a string array,
%   returns a string array by concatenating corresponding elements
%   of S1,S2, etc.  The inputs must all have the same size
%   (or any can be a scalar).

%   Copyright 2016 The MathWorks, Inc.

% Convert string arguments to cell arrays
args = varargin;
for idx = 1:numel(args)
    if isstring(args{idx})
        args{idx} = cellstr(args{idx});
    end
end

t = string(strcat(args{:}));
