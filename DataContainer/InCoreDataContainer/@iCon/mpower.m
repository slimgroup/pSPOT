function y = mpower(A,B,swp)
%^   Matrix power.
%
%   A^y is A to the y power.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.
   
%   http://www.cs.ubc.ca/labs/scl/spot

% unswap
if nargin == 3 && strcmp(swp,'swap')
    tmp = B;
    B = A;
    A = tmp;
    clear('tmp');
end

if size(A,1) ~= size(A,2)
   error('Operator must be square.');
end
if ~isscalar(B)
   error('Exponent must be a scalar.');
end

if isa(A,'iCon')
    y      = A;
    y.data = double(A) ^ double(B);
    
elseif isa(B,'iCon')
    y = A ^ double(B);
    
else
   error('Invalid parameters to mpower.');
end