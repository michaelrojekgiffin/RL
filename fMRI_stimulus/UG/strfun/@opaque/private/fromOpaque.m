function z = fromOpaque(x)

% Copyright 2006 The MathWorks, Inc.

z=x;

if isjava(z)
    z = char(z);
end

if isa(z,'opaque')
    error(message('MATLAB:fromOpaque:InvalidConversion', class( x )));
end