function y = rdivide(A,B,swp)
%./  Right array divide.
%   A./B denotes element-by-element division.  A and B
%   must have the same dimensions unless one is a scalar.
%   A scalar can be divided with anything.
%
%   C = RDIVIDE(A,B) is called for the syntax 'A ./ B' when A or B is an
%   object.
%
%   See also LDIVIDE

if nargin == 3 && strcmp(swp,'swap')
   temp = A;
   A = B;
   B = temp;
   clear temp;
end

% Case for both scalars
if isscalar(A) && isscalar(B)
    y = double(A) ./ double(B);
    return;
end

if isscalar(A)
    A = A*ones(size(B));
end
if isscalar(B)
    B = B*ones(size(A));
end

if ~isa(A,'iCon') % Right divide
    y      = B;
    y.data = double(A ./ double(B));
            
elseif ~isa(B,'iCon') % Left divide
    y      = A;
    y.data = double(double(A) ./ B);
    
else % Both data containers
    y      = A;
    y.data = double(A) ./ double(B);
    
    % Check for strict flag
    if A.strict || B.strict
       assert(all(A.imdims == B.imdims),...
           'Strict flag enforced. Implicit dimensions much match.')
    end
end