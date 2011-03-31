function result = transpose(x)
%.'   Operator tranpose.
%   A.' is the (non-conjugate) transpose of A.
%
%   transpose(A) is called for the syntax A.' when A is an operator.
%
%   See also opTranspose, opCTranspose, opSpot.ctranspose.
   
%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

% Check for ndims
assert(ndims(x) == 2, 'x must be 2D')

% Transpose
result        = x;
result.data   = transpose(x.data);
result.exdims = fliplr(x.exdims);
result.imdims = circshift(x.imdims,[0 ...
    DataContainer.collapsedDim(x.exdims,x.imdims) - 1]);
