function y = power(A,B)
%.^  Array power.
%   Z = X.^Y denotes element-by-element powers.  X and Y
%   must have the same dimensions unless one is a scalar. 
%   A scalar can operate into anything.

y = iCon(power(double(A),double(B)));

if isa(A,'dataContainer')
    y.imdims = A.imdims;
else
    y.imdims = B.imdims;
end