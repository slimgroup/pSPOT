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
    
    % Extract collapsed dimensions
    y.imdims = { size(A,1) B.imdims{2} };
            
elseif ~isa(B,'dataContainer') % Left multiply
    y        = A;
    y.data   = double(A)*B;
    y.exdims = size(y.data);
    
    % Extract collapsed dimensions
    y.imdims = { A.imdims{1} size(B,2) };
    
else % Both data containers
    y        = A;
    y.data   = double(A)*double(B);
    y.exdims = size(y.data);
    
    % Extract collapsed dimensions
    y.imdims = { A.imdims{1} B.imdims{2} };
end

end % mtimes