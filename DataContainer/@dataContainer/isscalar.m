function result = isscalar(x)
%ISSCALAR  True if operator is a scalar.
%
%   isscalar(x) returns true if x is a 1-by-1 operator.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

if (size(x,1) == 1) && (size(x,2) == 1)
   result = true;
else
   result = false;
end
