function y = minus(A,B,swp)
%+   Sum of two operators.
%
%   See also opSum, opSpot.minus.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

if nargin == 3 && strcmp(swp,'swap')
   temp = A;
   A = B;
   B = temp;
   clear temp;
end

% Case for both scalars
if isscalar(A) && isscalar(B)
    y = double(A) - double(B);
    return;
end

if isscalar(A)
    A = A*ones(size(B));
end
if isscalar(B)
    B = B*ones(size(A));
end

if ~isa(A,'iCon') % Right minus
    y      = B;
    y.data = double(A - double(B));
            
elseif ~isa(B,'iCon') % Left minus
    y      = A;
    y.data = double(double(A) - B);
    
else % Both data containers
    y      = A;
    y.data = double(A) - double(B);
    
    % Check for strict flag
    if A.strict || B.strict
       assert(all(A.imdims == B.imdims),...
           'Strict flag enforced. Implicit dimensions much match.')
    end
end