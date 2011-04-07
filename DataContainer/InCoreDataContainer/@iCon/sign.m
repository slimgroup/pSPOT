function y = sign(x)
%SIGN   Signum function.
%   For each element of X, SIGN(X) returns 1 if the element
%   is greater than zero, 0 if it equals zero and -1 if it is
%   less than zero.  For the nonzero elements of complex X,
%   SIGN(X) = X ./ ABS(X).
%
%   See also ABS.

y = x;
y.data = sign(double(x));