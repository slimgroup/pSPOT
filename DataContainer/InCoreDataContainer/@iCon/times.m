function y = times(A,B,swp)
%.*  Array multiply.
%   X.*Y denotes element-by-element multiplication.  X and Y
%   must have the same dimensions unless one is a scalar.
%   A scalar can be multiplied into anything.
%
%   C = TIMES(A,B) is called for the syntax 'A .* B' when A or B is an
%   object.
%
%   See also MTIMES.

if nargin == 3 && strcmp(swp,'swap')
   temp = A;
   A = B;
   B = temp;
   clear temp;
end

% Case for both scalars
if isscalar(A) && isscalar(B)
    y = double(A) .* double(B);
    return;
end

if isscalar(A)
    A = iCon(A*ones(size(B)));
end
if isscalar(B)
    B = iCon(B*ones(size(A)));
end

if ~isa(A,'iCon') % Right multiply
    y      = B;
    y.data = double(A .* double(B));
            
elseif ~isa(B,'iCon') % Left multiply
    y      = A;
    y.data = double(double(A) .* B);
    
else % Both data containers
    y      = A;
    y.data = double(A) .* double(B);
    
    % Check for strict flag
    if A.strict || B.strict
       assert(all(A.imdims == B.imdims),...
           'Strict flag enforced. Implicit dimensions much match.')
    end
end