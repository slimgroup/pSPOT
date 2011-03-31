function y = plus(A,B,swp)
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

if isscalar(A)
   A = dcInCore(A*ones(size(B)));
end
if isscalar(B)
   B = dcInCore(B*ones(size(A)));
end

if ~isa(A,'dataContainer') % Right multiply
    y = dcInCore(double(A + double(B)));
    
    y.imdims = B.imdims;
            
elseif ~isa(B,'dataContainer') % Left multiply
    y = dcInCore(double(double(A) + B));
    
    y.imdims = A.imdims;
    
else % Both data containers
    y = dcInCore(double(A) + double(B));
    
    % Check for strict flag
    if A.strict || B.strict
       assert(all(A.imdims == B.imdims),...
           'Strict flag enforced. Implicit dimensions much match.')
    end
    y.imdims = A.imdims;
end