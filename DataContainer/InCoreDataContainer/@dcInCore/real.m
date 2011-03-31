function y = real(x)
%REAL  Complex real part.
%
%   opReal(A) returns an operator comprised of the real part of A.
%
%   See also opReal, opSpot.imag.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

   y = dcInCore(real(double(x)));
   y.imdims = x.imdims;