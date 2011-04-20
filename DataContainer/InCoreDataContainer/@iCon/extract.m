function y = extract(x)
%EXTRACT    Extract metadata of data container
%
%   extract(x) returns an empty data container x containing all the vital
%   metadata stores in x minus the actual data itself, as well as having
%   its explicit dimensions set to zero.
%
%   See also: inject

y        = x;
y.data   = [];
y.exdims = [0 0];