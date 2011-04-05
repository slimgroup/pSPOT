function y = ldivide(A,B)
%.\  Left array divide.
%   A.\B denotes element-by-element division.  A and B
%   must have the same dimensions unless one is a scalar.
%   A scalar can be divided with anything.
%
%   C = LDIVIDE(A,B) is called for the syntax 'A .\ B' when A or B is an
%   object.
%
%   See also RDIVIDE

y = iCon(ldivide(double(A),double(B)));

if isa(A,'dataContainer')
    y.imdims = A.imdims;
else
    y.imdims = B.imdims;
end