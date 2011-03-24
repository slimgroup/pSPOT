function y = ndims(x)
%   N = NDIMS(X) returns the number of dimensions in the array X.
%   The number of dimensions in an array is always greater than
%   or equal to 2.  Trailing singleton dimensions are ignored.
%   Put simply, it is LENGTH(SIZE(X)).

y = length(size(x));