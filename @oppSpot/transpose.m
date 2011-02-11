function result = transpose(A)
%.'   Operator tranpose for pSpot
%   (the sole purpose of this operator is to direct transposed 
%   pSpot operators to the correct mtimes)
%
%   A.' is the (non-conjugate) transpose of A.
%
%   transpose(A) is called for the syntax A.' when A is an operator.
%
%   See also opTranspose, opCTranspose, opSpot.ctranspose.
   
%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

result = oppTranspose(A);
