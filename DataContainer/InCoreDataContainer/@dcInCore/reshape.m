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

% Check for the collapsibility of reshape
% Do the calculation
imdims  = x.imdims;
exdims  = [varargin{:}];
result1 = -1;
collapsed_dims = 1;
for i = 1:length(imdims)
    collapsed_dims = collapsed_dims * imdims(i);
    if  collapsed_dims == exdims(1)
        result1 = i + 1;
        break;
    end
end

exdims  = x.imdims;
imdims  = [varargin{:}];
result2 = -1;
collapsed_dims = 1;
for i = 1:length(imdims)
    collapsed_dims = collapsed_dims * imdims(i);
    if  collapsed_dims == exdims(1)
        result2 = i + 1;
        break;
    end
end

% Boolean
result = result1 > 0 || result2 > 0;

assert(result, ...
    'Reshape dimensions must be collapesed or multiples of implicit dimension')

% Reshape
y = dcInCore(reshape(x.data,varargin{:}));
y.imdims = x.imdims;