function tf = contains(s, pattern, varargin)
%CONTAINS True if pattern is found in text.
%   TF = CONTAINS(S,PATTERN) returns true if CONTAINS finds PATTERN in any
%   element of string array S. TF is the same size as S.
%
%   S can be a string array, a character vector, or a cell array of
%   character vectors. So can PATTERN. PATTERN and S need not be the same
%   size. If PATTERN is nonscalar, CONTAINS returns true if it finds any
%   element of PATTERN in S.
% 
%   TF = CONTAINS(S,PATTERN,'IgnoreCase',IGNORE) ignores case when searching 
%   for PATTERN in S if IGNORE is true. The default value of IGNORE is false.
% 
%   Examples
%       S = string('data.tar.gz');
%       P = string('tar');
%       contains(S,P)                   returns  1
%
%       S = string({'abstracts.docx','data.tar.gz'});
%       P = 'tar';         
%       contains(S,P)                   returns  [0 1]
%
%       S = string('data.tar.gz');
%       P = {'docx','tar'};
%       contains(S,P)                   returns  1
%
%       S = string({'DATA.TAR.GZ','SUMMARY.PPT'});
%       P = string('tar');
%       contains(S,P,'IgnoreCase',true) returns  [1 0]
%
%   See also endsWith, startsWith.

%   Copyright 2015-2016 The MathWorks, Inc.

    narginchk(2, inf);

    if ~ischar(s) && ~iscellstr(s) && ~isstring(s)
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    end

    try
        stringS = string(s);
        tf = stringS.contains(pattern, varargin{:});
    catch E
        throw(E)
    end
end
