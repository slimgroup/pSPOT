function y = mldivide(A,B,swp)
%\   Backslash or left matrix divide.
%   A\B is the matrix division of A into B, which is roughly the
%   same as INV(A)*B , except it is computed in a different way.
%   If A is an N-by-N matrix and B is a column vector with N
%   components, or a matrix with several such columns, then
%   X = A\B is the solution to the equation A*X = B. A warning
%   message is printed if A is badly scaled or nearly singular.
%   A\EYE(SIZE(A)) produces the inverse of A.
%
%   If A is an M-by-N matrix with M < or > N and B is a column
%   vector with M components, or a matrix with several such columns,
%   then X = A\B is the solution in the least squares sense to the
%   under- or overdetermined system of equations A*X = B. The
%   effective rank, K, of A is determined from the QR decomposition
%   with pivoting. A solution X is computed which has at most K
%   nonzero components per column. If K < N this will usually not
%   be the same solution as PINV(A)*B.  A\EYE(SIZE(A)) produces a
%   generalized inverse of A.

% unswap
if nargin == 3 && strcmp(swp,'swap')
    tmp = B;
    B = A;
    A = tmp;
    clear('tmp');
end

if isscalar(A)
    y = A .\ B;
    
elseif ~isa(A,'iCon')
    y        = B;
    y.data   = double( A \ double(B) );
    y.exdims = size(y.data);
    
    % Extract collapsed dimensions & permutation
    y.imdims = { size(A,2) B.imdims{2} };
    y.perm   = B.perm;
    
    % Check for spot ms and ns
    if isa(A,'opSpot')
        y.imdims{1} = A.ns;
    end
    
elseif ~isa(B,'iCon')
    y        = A;
    y.data   = double( double(A) \ B );
    y.exdims = size(y.data);
    
    % Extract collapsed dimensions & permutation
    y.imdims = { A.imdims{2} size(B,2) };
    y.perm   = A.perm;
    
    % Check for spot ms and ns
    if isa(A,'opSpot')
        y.imdims{2} = B.ns;
    end
    
else % Both data containers
    y        = A;
    y.data   = double(A) \ double(B);
    y.exdims = size(y.data);
    
    % Extract collapsed dimensions
    y.imdims = { A.imdims{2} B.imdims{2} };
    y.perm   = A.perm;
end

end % mldivide