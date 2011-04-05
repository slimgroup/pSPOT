function y = power(A,B)
%.^  Array power.
%   Z = X.^Y denotes element-by-element powers.  X and Y
%   must have the same dimensions unless one is a scalar. 
%   A scalar can operate into anything.
%
%   C = POWER(A,B) is called for the syntax 'A .^ B' when A or B is an
%   object.
%
%   See also MPOWER

y = iCon(power(double(A),double(B)));

if isa(A,'dataContainer')
    y.imdims = A.imdims;
else
    y.imdims = B.imdims;
end