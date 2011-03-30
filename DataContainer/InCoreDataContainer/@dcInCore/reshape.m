function y = reshape(x,varargin)
%RESHAPE    Reshape data container object to desired shape
%
%   reshape(X,N1,N2,...,N) reshapes data container X into the
%   dimensions defined as [N1,N2,...,N]. Note that the number of elements
%   must be conserved.
%
%   Always keep in mind that reshape on distributed arrays always conserve
%   the elements locally on the labs, ie. There will be no communication
%   between labs. Therefore, the local parts size after reshaping has to be
%   the same locally. This is generally not an issue if you preserve the
%   size of the distributed dimension. Or some special symmetrical
%   distribution scheme is used.
%
%   See also: unvec, vec, double

y = dcInCore(reshape(x.data,varargin{:}));
y.imdims = x.imdims;