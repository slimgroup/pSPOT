function y = imag(x)
%IMAG  Complex imaginary part.
%
%   imag(A) returns a new operator comprised of imaginary part of A.
%
%   See also opSpot.real, opImag.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

y = iCon(imag(double(x)));
y.imdims = x.imdims;