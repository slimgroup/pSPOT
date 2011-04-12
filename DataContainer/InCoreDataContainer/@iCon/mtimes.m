function y = mtimes(A,B,swp)
% unswap
if nargin == 3 && strcmp(swp,'swap')
    tmp = B;
    B = A;
    A = tmp;
    clear('tmp');
end

% Multiply
if ~isa(A,'dataContainer') % Right multiply
    y        = B;
    y.data   = A*double(B);
    y.exdims = size(y.data);
    
    % Extract collapsed dimensions & permutation
    y.imdims = { size(A,1) B.imdims{2} };
    y.perm   = B.perm;
    
    % Check for spot ms and ns
    if isa(A,'opSpot')
        y.imdims{1} = A.ms;
    end
    
elseif ~isa(B,'dataContainer') % Left multiply
    y        = A;
    y.data   = double(A)*B;
    y.exdims = size(y.data);
    
    % Extract collapsed dimensions & permutation
    y.imdims = { A.imdims{1} size(B,2) };
    y.perm   = A.perm;
    
    % Check for spot ms and ns
    if isa(A,'opSpot')
        y.imdims{2} = A.ns;
    end
    
else % Both data containers
    y        = A;
    y.data   = double(A)*double(B);
    y.exdims = size(y.data);
    
    % Extract collapsed dimensions
    y.imdims = { A.imdims{1} B.imdims{2} };
end

end % mtimes