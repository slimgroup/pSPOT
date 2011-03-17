function x = univec(x)
%UNIVEC     Removes implicit vectorization of data container x
%
%   univec(x) sets the implicit vectorization status of x to false.
%
%   See also: ivec, vec, unvec, double
x.ivec = false;