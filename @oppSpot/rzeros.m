function y = rzeros(A)
%RZEROS Distributed zero vector in operator range
%
%   y = rzeros(A) generates a zero vector with the size of the operator
%   range, distributed accordingly to the operators needs so that A'*y is
%   a valid operation.

