function y = power(A,B)
%.^  Array power.
%   Z = X.^Y denotes element-by-element powers.  X and Y
%   must have the same dimensions unless one is a scalar.
%   A scalar can operate into anything.
%
%   See also MPOWER

y = piCon(power(double(A),double(B)));

if isa(A,'dataContainer')
    y.imdims    = A.imdims;
    y.imcoddims = A.imcoddims;
    y.imcodpart = A.imcodpart;
else
    y.imdims    = B.imdims;
    y.imcoddims = B.imcoddims;
    y.imcodpart = B.imcodpart;
end