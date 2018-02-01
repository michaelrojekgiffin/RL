function i = strmatch(str,strs,flag)
%STRMATCH Find possible matches for character arrays.
%   IDX = STRMATCH(STR, STRARRAY) looks through the rows of the character
%   array or cell array of character vectors STRARRAY to find character
%   vector that begin with STR, and returns the matching row indices. Any
%   trailing space characters in STR or STRARRAY are ignored when matching.
%   STRMATCH is fastest when STRARRAY is a character array.
%
%   IDX = STRMATCH(STR, STRARRAY, 'exact') compares STR with each row of
%   STRARRAY, looking for an exact match of the entire character vector.
%   Any trailing space characters in STR or STRARRAY are ignored when
%   matching.
%
%   Examples
%     idx = strmatch('max',strvcat('max','minimax','maximum'))
%   returns idx = [1; 3] since rows 1 and 3 begin with 'max', and
%     idx = strmatch('max',strvcat('max','minimax','maximum'),'exact')
%   returns idx = 1, since only row 1 matches 'max' exactly.
%   
%   STRMATCH will be removed in a future release. Use startsWith or
%   validatestring instead.
%
%   See also startsWith, VALIDATESTRING, CONTAINS, STRFIND, STRCMP,
%   STRNCMP, REGEXP

%   Copyright 1984-2016 The MathWorks, Inc.

% The cell array implementation is in @cell/strmatch.m

narginchk(2,3);

[m,n] = size(strs);
len = numel(str);

if (nargin==3)
    exactMatch = true;
    if ~ischar(flag)
        warning(message('MATLAB:strmatch:InvalidFlagType'));
    elseif ~strcmpi(flag,'exact')
        warning(message('MATLAB:strmatch:InvalidFlag', flag, flag));
    end
else
    exactMatch = false;
end

% Special treatment for empty STR or STRS to avoid
% warnings and error below
if len==0
    str = reshape(str,1,len);
end 
if n==0
    strs = reshape(strs,max(m,1),n);
    [m,n] = size(strs);
end

if len > n
    i = [];
else
    if exactMatch && len < n % if 'exact' flag, pad str with blanks or nulls
        [~,strn] = size(str);
        if strn ~= len
            error(message('MATLAB:strmatch:InvalidShape'));
        else
            % Use nulls if anything in the last column is a null.
            null = char(0); 
            space = ' ';
            if ~isempty(strs) && any(strs(:,end)==null), 
                str = [str null(ones(1,n-len))];
            else
                str = [str space(ones(1,n-len))];
            end
            len = n;
        end
    end

    mask = true(m,1); 
    % walk from end of strs array and search for row starting with str.
    for outer = 1:m
        for inner = 1:len
            if (strs(outer,inner) ~= str(inner))
                mask(outer) = false;
                break; % exit matching this row in strs with str.
            end   
        end
    end 
    i = find(mask);
end
