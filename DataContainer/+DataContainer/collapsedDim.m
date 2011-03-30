function y = collapsedDim(x)
%COLLAPSEDDIM   Search Algorithm for the Fundamental Theorem of Henryk
%
%   collapsedDim(exdims,imdims) returns the implicit dimension that is the
%   contiguous dimension in the second dimension of the explicit dimensions
%   ie. If x has implicit dimensions n1 x n2 x n3 and explicit dimensions
%   n1 x n2*n3, then collapsedDim will return 2, indicating that n2 is the
%   contiguous dimension of the collapsed dimensions n2*n3
%
%   Note: This algorithm only applies if explicit dimensions are 2D

% Check nargin
assert(nargin == 1, 'Must have exactly 1 input arguments')

% Setup variables
exdims = x.exdims;
imdims = x.imdims;

% Check 2D-ness of exdims
assert(length(exdims) == 2, 'Exdims must be 2D')

% Do the calculation
y = 2;
collapsed_dims = 1;
for i = 1:length(imdims)
    collapsed_dims = collapsed_dims * imdims(i);
    if  collapsed_dims == exdims(1)
        y = i + 1;
        break;
    end
end