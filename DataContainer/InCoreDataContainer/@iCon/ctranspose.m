function result = ctranspose(x)
%'  Complex conjugate tranpose.
%   A' is the complex conjugate transpose of A.
%
%   ctranspose(A) is called for the syntax A' when A is an operator.
%
%   See also opCTranspose, opTranspose, opSpot.transpose.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

% Conjugate Transpose
result        = x;
result.data   = ctranspose(double(x));
result.exdims = fliplr(x.exdims);
result.imdims = fliplr(x.imdims);
result.perm   = fliplr(x.perm);