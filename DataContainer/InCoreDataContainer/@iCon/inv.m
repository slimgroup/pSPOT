function y = inv(x)
%INV    Matrix inverse.
%   INV(X) is the inverse of the square matrix X.
%   A warning message is printed if X is badly scaled or
%   nearly singular.

y        = x;
y.data   = inv(double(x));
y.exdims = fliplr(x.exdims);
y.imdims = fliplr(x.imdims);
y.perm   = fliplr(x.perm);