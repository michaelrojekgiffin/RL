function len = strlength(s)
%STRLENGTH Lengths of text elements.
%   L = STRLENGTH(STR) returns a numeric array where each element is the
%   number of characters in the corresponding string in STR. L is the same
%   size as STR.
% 
%   Example:
%       STR = string('data.xlsx');
%       strlength(STR)      
%
%       returns  
%
%            9
% 
%   Example:
%       STR = string({'funds.xlsx';'demo.ppt'});
%       strlength(STR)      
%
%       returns  
%
%           10
%            8
%
%   See also LENGTH.

%   Copyright 2014-2016 The MathWorks, Inc.

    narginchk(1, 1);

    try
        if iscell(s)
            len = zeros(size(s));
            for idx = 1:numel(s)
                len(idx) = charlength(s{idx});
            end
        else
            len = charlength(s);
        end
    catch E
        throw(E)
    end
end

function clen = charlength(s)

    if ~ischar(s) || (~isempty(s) && ~isrow(s))
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    end
    clen = numel(s);

end
