function y = mpower(A,B,swp)
%^   Matrix power.
%
%   A^y is A to the y power.
%
%   Note: A must be square, and y must be scalar

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