function y = double(x)
%DOUBLE     Returns the data contained in the data container
%
%   double(A) returns the underlying data contained in the data container
%   If x is implicitly vectorized, double will make it explicit.
%
%   See also: vec, reshape

y = x.data;