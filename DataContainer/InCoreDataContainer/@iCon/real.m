function y = real(x)
%REAL  Complex real part.
%
%   real(A) returns an operator comprised of the real part of A.
%
%   See also imag.

y      = x;
y.data = real(double(x));